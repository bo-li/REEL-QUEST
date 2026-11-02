#!/usr/bin/env python3
"""
REEL-QUEST (Rare-Earth Extraction Ligands - Quantitative Understanding of Elemental Selectivity Trends)
Interactive Dashboard - Streamlit-based web interface for browsing and visualizing ligand extraction data
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import numpy as np
from scipy import stats
import sys
import os
import re
import hashlib

# Add parent directory to path for rel_quest import
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from rel_quest import LigandDatabase
from rel_quest.schema import TREND_DISPLAY_ORDER
from rel_quest.references import doi_url, primary_doi
# Shared constants and chemistry/data helpers live in rel_quest.plotting —
# the single source of truth also used by non-Streamlit scripts.
from rel_quest.plotting import (
    LANTHANIDES,
    ATOMIC_NUMBERS,
    EXTRACTANT_TYPE_COLORS,
    CHEMICAL_FAMILY_COLORS,
    SOLVENT_COLORS,
    SOLVENT_CATEGORY_COLORS,
    LIGAND_ACIDITY_COLORS,
    FONT_SIZES_NORMAL,
    FONT_SIZES_FULLSCREEN,
    load_solvent_categories,
    get_solvent_category,
    classify_ligand_acidity,
    extract_element_data,
    extract_separation_factor_data,
    split_bounds,
    category_order,
    bound_text,
    bound_trace,
    BOUND_COLOR,
)
from rel_quest.uncertainty import is_censored


# LLM agent (optional — graceful fallback if requests not installed)
try:
    from rel_quest.agent import RELQuestAgentStreaming, check_ollama_connection, list_available_models, DEFAULT_MODEL
    from rel_quest.llm_backends import (
        OllamaBackend,
        OpenAICompatibleBackend,
        OLLAMA_BASE_URL,
        OPENROUTER_BASE_URL,
        get_config,
        is_free_model,
    )
    AGENT_AVAILABLE = True
except ImportError:
    AGENT_AVAILABLE = False

# Try to import RDKit for molecular structures. Only the structure panel needs
# it, so a failure is reported there rather than as a banner on every page —
# and never before st.set_page_config, which must come first.
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

# Page configuration
st.set_page_config(
    page_title="REEL-QUEST - Rare-Earth Extraction Ligand Database",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Lanthanide element data for periodic table display
LANTHANIDE_INFO = {
    'La': {'atomic_number': 57, 'name': 'Lanthanum', 'color': '#B8E6F0'},
    'Ce': {'atomic_number': 58, 'name': 'Cerium', 'color': '#B8E6F0'},
    'Pr': {'atomic_number': 59, 'name': 'Praseodymium', 'color': '#B8E6F0'},
    'Nd': {'atomic_number': 60, 'name': 'Neodymium', 'color': '#B8E6F0'},
    'Pm': {'atomic_number': 61, 'name': 'Promethium', 'color': '#B8E6F0'},
    'Sm': {'atomic_number': 62, 'name': 'Samarium', 'color': '#B8E6F0'},
    'Eu': {'atomic_number': 63, 'name': 'Europium', 'color': '#B8E6F0'},
    'Gd': {'atomic_number': 64, 'name': 'Gadolinium', 'color': '#B8E6F0'},
    'Tb': {'atomic_number': 65, 'name': 'Terbium', 'color': '#B8E6F0'},
    'Dy': {'atomic_number': 66, 'name': 'Dysprosium', 'color': '#B8E6F0'},
    'Ho': {'atomic_number': 67, 'name': 'Holmium', 'color': '#B8E6F0'},
    'Er': {'atomic_number': 68, 'name': 'Erbium', 'color': '#B8E6F0'},
    'Tm': {'atomic_number': 69, 'name': 'Thulium', 'color': '#B8E6F0'},
    'Yb': {'atomic_number': 70, 'name': 'Ytterbium', 'color': '#B8E6F0'},
    'Lu': {'atomic_number': 71, 'name': 'Lutetium', 'color': '#B8E6F0'}
}

# DOI of the paper that measured a ligand's values (Ligand Browser column).
# Since schema v0.3 `reference.DOI` is already the bare primary DOI; the shared
# helper also reads older generations, where the field could hold a
# comma-joined pair ending in SAFE's SAFE-MolGen curation tag
# (10.1021/jacs.5c19738 — a curation tag, not a measurement source).
def _get_primary_doi(ligand):
    """First measuring-paper DOI across the ligand's conditions, or None."""
    for cond in (getattr(ligand, 'experimental_conditions', None) or []):
        if isinstance(cond, dict):
            doi = primary_doi(cond.get('reference'))
            if doi:
                return doi
    return None


#=============================================================================
# Helper Functions
#=============================================================================

@st.cache_resource
def load_database():
    """Load the ligand database (cached for performance).

    Path resolution is delegated to LigandDatabase, which resolves the current
    canonical file relative to the rel_quest package (CWD-independent) and
    warns if it has to fall back to an older generation. RELQUEST_DB names
    another file to preview it (e.g. an unpromoted candidate); the page then
    says which file it shows.
    """
    return LigandDatabase(get_config("RELQUEST_DB") or None)


def render_molecule(smiles, size=(400, 300)):
    """Render 2D molecular structure from SMILES using RDKit"""
    if not RDKIT_AVAILABLE:
        return None

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, size=size)
        return img
    except Exception as e:
        st.warning(f"Error rendering molecule: {e}")
        return None


def get_font_sizes():
    """
    Get font sizes based on fullscreen_fonts mode.

    Returns:
        dict: Font size configuration (FONT_SIZES_NORMAL or FONT_SIZES_FULLSCREEN)
    """
    if st.session_state.get('fullscreen_fonts', False):
        return FONT_SIZES_FULLSCREEN
    else:
        return FONT_SIZES_NORMAL


def render_element_selector(state_key, button_prefix, label, caption_prefix):
    """
    Render the periodic-table lanthanide selector: 15 buttons, La through Lu.

    The selected element is highlighted (type="primary"); clicking any other
    element writes it to session state and reruns.

    Args:
        state_key: st.session_state key holding the current selection
        button_prefix: prefix for widget keys, must be unique per selector
        label: markdown heading rendered above the buttons
        caption_prefix: prefix for the caption below (e.g. "Selected: ")

    Returns:
        str: the currently selected element symbol
    """
    st.markdown(label)

    cols = st.columns(15)
    for idx, element in enumerate(LANTHANIDES):
        with cols[idx]:
            atomic_num = LANTHANIDE_INFO[element]['atomic_number']
            button_label = f"{atomic_num}\n{element}"
            help_text = f"{LANTHANIDE_INFO[element]['name']} (Z={atomic_num})"
            key = f"{button_prefix}_{element}"

            if st.session_state[state_key] == element:
                st.button(button_label, key=key, type="primary",
                          help=help_text, use_container_width=True)
            elif st.button(button_label, key=key, help=help_text,
                           use_container_width=True):
                st.session_state[state_key] = element
                st.rerun()

    selected = st.session_state[state_key]
    st.caption(f"{caption_prefix}{LANTHANIDE_INFO[selected]['name']} ({selected})")
    return selected


def _empty_plot_figure(message, title=None, hide_axes=False):
    """
    Build the "no data" placeholder figure shared by the distribution plots.

    Args:
        message: centred text to display
        title: optional figure title; when given, height/template are set too
        hide_axes: hide tick labels and gridlines (used when there are no axes
            worth showing at all)

    Returns:
        plotly.graph_objects.Figure
    """
    fig = go.Figure()
    font_sizes = get_font_sizes()
    fig.add_annotation(
        text=message,
        xref="paper", yref="paper",
        x=0.5, y=0.5, showarrow=False,
        font=dict(size=font_sizes['annotation'], color="gray")
    )
    if title is not None:
        layout = dict(title=title, height=500, template="plotly_white")
        if hide_axes:
            layout['xaxis'] = dict(showticklabels=False, showgrid=False)
            layout['yaxis'] = dict(showticklabels=False, showgrid=False)
        fig.update_layout(**layout)
    return fig


def _category_color_map(x_axis_option):
    """
    Return the colour palette for a categorical x-axis grouping.

    Callers supply their own fallback colour via .get(category, default) —
    the two distribution plots deliberately use different defaults.
    """
    if x_axis_option == "Ligand Types":
        # Chemical_Family palette wins; legacy extractant-type colours fill gaps.
        return {**EXTRACTANT_TYPE_COLORS, **CHEMICAL_FAMILY_COLORS}
    elif x_axis_option == "Solvent Categories":
        return SOLVENT_CATEGORY_COLORS
    elif x_axis_option == "Ligand Acidity":
        return LIGAND_ACIDITY_COLORS
    else:
        return SOLVENT_COLORS


def calculate_condition_r2(distribution_ratios):
    """Calculate R² for a single experimental condition"""
    elements = []
    z_vals = []
    log_d_vals = []

    for elem in LANTHANIDES:
        if elem in distribution_ratios and not is_censored(distribution_ratios[elem]):   # a bound is no point
            d_val = distribution_ratios[elem].get('value', 0)
            if d_val > 0:
                elements.append(elem)
                z_vals.append(ATOMIC_NUMBERS[elem])
                log_d_vals.append(np.log10(d_val))

    if len(z_vals) < 4:
        return None, None

    slope, intercept, r_value, p_value, std_err = stats.linregress(z_vals, log_d_vals)
    r_squared = r_value ** 2

    return r_squared, slope


def condition_bounds_text(distribution_ratios):
    """'Nd > 316; La < 0.01': the bounds of one condition ('' when it has none)."""
    return '; '.join(f"{el} {d['censor']} {10 ** d['log_value']:.3g}"
                     for el, d in distribution_ratios.items() if is_censored(d))


def co_extractants_text(condition):
    """'TODGA 0.005 M': the other extractants of a mixed-solvent condition ('' when it has none)."""
    return '; '.join(f"{x['Name']}" + (f" {x['Concentration_M']:g} M" if x.get('Concentration_M') is not None else '')
                     for x in condition.get('Co_Extractants') or [])


def has_value_in_range(ligand, logd_min, logd_max):
    """Whether the ligand has a measured log D in [logd_min, logd_max]; a bound is not a value."""
    return any(not is_censored(data) and data.get('log_value') is not None
               and logd_min <= data['log_value'] <= logd_max
               for cond in ligand.experimental_conditions
               for data in cond.get('distribution_ratios', {}).values())


def _element_bound_hover(row, element):
    return (f"<b>{row['ligand_name']}</b><br>{bound_text(f'D[{element}]', row['censor'], row['log_d'])}<br>"
            f"Solvent: {row['solvent']}<br>Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"
            f"Extractant Conc: {row['extractant_conc']:.3f} M<br>Temperature: {row['temperature']:.1f} °C")


def _sf_bound_hover(row, element1, element2):
    def log_d(value, censor):
        # censor is '>', '<', or missing (None, or NaN once pandas has stored the column)
        return f"{censor + ' ' if isinstance(censor, str) else ''}{value:.2f}"
    text = (f"<b>{row['ligand_name']}</b><br>{bound_text(f'SF[{element1}/{element2}]', row['censor'], row['log_sf'])}<br>"
            f"logD[{element1}]: {log_d(row['log_d_elem1'], row['censor_elem1'])}<br>"
            f"logD[{element2}]: {log_d(row['log_d_elem2'], row['censor_elem2'])}<br>"
            f"Solvent: {row['solvent']}<br>Acid: {row['acid_name']}")
    if row['acid_conc'] is not None and not pd.isna(row['acid_conc']):
        text += f" ({row['acid_conc']:.2f} M)"
    return text


def plot_selectivity_trend(ligand, condition_id=None):
    """Create interactive selectivity trend plot (log(D) vs atomic number)"""
    # Get data for specific condition or aggregated
    if condition_id:
        conditions = [c for c in ligand.experimental_conditions if c['condition_id'] == condition_id]
        if not conditions:
            return None
        cond = conditions[0]
        d_ratios = cond.get('distribution_ratios', {})
        title_suffix = f"<br><sub>Condition: {condition_id[:60]}...</sub>"
    else:
        # Aggregate across all conditions (median)
        # Note: get_aggregated_d_values already returns dicts with 'value' keys
        d_ratios = ligand.get_aggregated_d_values(aggregation='median', lanthanides_only=True)
        title_suffix = "<br><sub>Aggregated (median across all conditions)</sub>"

    # Extract data points
    elements = []
    z_vals = []
    d_vals = []
    log_d_vals = []
    sigma_log_d_vals = []
    bounds = []   # (Z, log bound, censor, element): drawn, never fitted

    for elem in LANTHANIDES:
        if elem in d_ratios:
            d_entry = d_ratios[elem]
            if is_censored(d_entry):
                bounds.append((ATOMIC_NUMBERS[elem], d_entry['log_value'], d_entry['censor'], elem))
                continue
            # Handle both dict format {'value': x} and plain numeric values
            if isinstance(d_entry, dict):
                d_val = d_entry.get('value', 0)
                sigma = d_entry.get('sigma_log_D', 0.033)
            else:
                d_val = d_entry
                sigma = 0.033

            # Ensure d_val is numeric and positive
            if isinstance(d_val, (int, float)) and d_val > 0:
                elements.append(elem)
                z_vals.append(ATOMIC_NUMBERS[elem])
                d_vals.append(d_val)
                log_d_vals.append(np.log10(d_val))
                sigma_log_d_vals.append(sigma if sigma else 0.033)

    if len(z_vals) < 2:
        return None

    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(z_vals, log_d_vals)
    r_squared = r_value ** 2

    # Determine trend label and marker color
    sel = getattr(ligand, 'selectivity', None) or {}
    stored_trend = sel.get('trend_classification')
    if stored_trend in ('Peaked', 'U_Shaped', 'Ascending', 'Reverse', 'Flat', 'Mixed'):
        trend = stored_trend
    elif r_squared < 0.6:
        trend = 'Flat/Mixed'
    elif slope > 0.02:
        trend = 'Ascending'
    elif slope < -0.02:
        trend = 'Reverse'
    else:
        trend = 'Flat'

    trend_colors = {
        'Ascending': 'green', 'Reverse': 'red',
        'Peaked': '#E76F51', 'U_Shaped': '#457B9D',
    }
    marker_color = trend_colors.get(trend, 'gray')

    fig = go.Figure()
    z_range = np.linspace(min(z_vals) - 0.5, max(z_vals) + 0.5, 60)

    # Both the trend line and the CI band pivot on the centroid of the
    # plotted points. The stored slope is the mean of per-condition fits
    # and the stored CI is its bootstrap interval over conditions
    # (rel_quest.uncertainty.bootstrap_trend_ci), so in the aggregated view
    # draw that slope and that band. In a single-condition view the stored
    # CI does not describe this condition, so show only its own OLS fit.
    # (The band used to be drawn as ci * Z + intercept, pivoting at Z = 0,
    # which made its height at Z ~ 64 equal to 64 x the CI width — several
    # dex for most ligands.)
    z_bar = float(np.mean(z_vals))
    y_bar = float(np.mean(log_d_vals))
    stored_slope = sel.get('slope')
    ci_lower = sel.get('slope_CI_lower')
    ci_upper = sel.get('slope_CI_upper')
    aggregated = condition_id is None

    if aggregated and ci_lower is not None and ci_upper is not None:
        y_lower = y_bar + ci_lower * (z_range - z_bar)
        y_upper = y_bar + ci_upper * (z_range - z_bar)
        fig.add_trace(go.Scatter(
            x=np.concatenate([z_range, z_range[::-1]]).tolist(),
            y=np.concatenate([y_upper, y_lower[::-1]]).tolist(),
            fill='toself',
            fillcolor='rgba(31,119,180,0.15)',
            line=dict(color='rgba(0,0,0,0)'),
            name='95% CI of trend slope (bootstrap over conditions)',
            showlegend=True,
            hoverinfo='skip',
        ))

    if aggregated and stored_slope is not None:
        line_slope, line_name = stored_slope, 'Trend slope (mean of per-condition fits)'
    else:
        # OLS through the plotted points; passes through the centroid, so
        # this equals slope * z + intercept.
        line_slope, line_name = slope, 'Linear fit'
    fig.add_trace(go.Scatter(
        x=z_range.tolist(),
        y=(y_bar + line_slope * (z_range - z_bar)).tolist(),
        mode='lines',
        line=dict(color=marker_color, width=2, dash='dash'),
        name=line_name,
        hoverinfo='skip',
    ))

    # Peak element marker for Peaked/U_Shaped
    peak_elem = sel.get('peak_element')
    if trend in ('Peaked', 'U_Shaped') and peak_elem and peak_elem in z_vals:
        # Highlight the peak/trough data point with a diamond marker
        peak_idx = elements.index(peak_elem)
        label = 'Peak' if trend == 'Peaked' else 'Trough'
        fig.add_trace(go.Scatter(
            x=[z_vals[peak_idx]],
            y=[log_d_vals[peak_idx]],
            mode='markers+text',
            marker=dict(size=16, color=marker_color, symbol='diamond',
                        line=dict(width=2, color='white')),
            text=[peak_elem],
            textposition='top center',
            name=f'{label} ({peak_elem})',
            hovertemplate=f'<b>{label}: {peak_elem}</b><br>log(D) = {log_d_vals[peak_idx]:.3f}<extra></extra>',
        ))

    # Data points with error bars
    fig.add_trace(go.Scatter(
        x=z_vals,
        y=log_d_vals,
        mode='markers',
        marker=dict(size=12, color=marker_color, line=dict(width=1, color='white')),
        error_y=dict(type='data', array=sigma_log_d_vals, visible=True,
                     color=marker_color, thickness=1.5, width=4),
        text=elements,
        customdata=d_vals,
        hovertemplate='<b>%{text}</b><br>Atomic #: %{x}<br>D = %{customdata:.3f}<br>log(D) = %{y:.3f}<extra></extra>',
        name='Data'
    ))
    if bounds:
        fig.add_trace(bound_trace([b[0] for b in bounds], [b[1] for b in bounds], [b[2] for b in bounds],
                                  [f"<b>{b[3]}</b><br>{bound_text('D', b[2], b[1])}<br>not in the fit"
                                   for b in bounds], marker_color, showlegend=True))

    fig.update_layout(
        title=f"Selectivity Trend: {ligand.name} ({trend}){title_suffix}",
        xaxis_title="Atomic Number",
        yaxis_title="log(D)",
        xaxis=dict(dtick=1, range=[56.5, 71.5]),
        hovermode='closest',
        showlegend=True,
        height=450,
        template="plotly_white"
    )

    return fig


def plot_ligand_type_pie(ligands):
    """Horizontal bar of Subtype counts (21 categories, sorted by count desc).

    Kept under the legacy name so existing callers don't have to change.
    A bar chart scales to 21 categories far better than a pie — every label
    stays readable and exact counts are visible without hovering.
    """
    type_counts = {}
    for lig in ligands:
        type_name = getattr(lig, 'chemical_family', None) or lig.extractant_type or 'Unknown'
        type_counts[type_name] = type_counts.get(type_name, 0) + 1

    # Sort ascending so the largest bar appears at the TOP (horizontal-bar idiom).
    sorted_types = sorted(type_counts.items(), key=lambda x: x[1])
    labels = [t[0] for t in sorted_types]
    values = [t[1] for t in sorted_types]
    total = sum(values)

    colors = [CHEMICAL_FAMILY_COLORS.get(label, '#999999') for label in labels]
    percents = [v / total * 100 if total else 0 for v in values]
    text_labels = [f"{v}  ({p:.1f}%)" for v, p in zip(values, percents)]

    fig = go.Figure(data=[go.Bar(
        x=values,
        y=labels,
        orientation='h',
        marker=dict(color=colors, line=dict(color='white', width=1)),
        text=text_labels,
        textposition='outside',
        textfont=dict(size=11),
        cliponaxis=False,
        hovertemplate='<b>%{y}</b><br>Count: %{x}<br>Percentage: %{customdata:.1f}%<extra></extra>',
        customdata=percents,
    )])

    fig.update_layout(
        title={
            'text': f"Ligand Subtype Distribution ({total} total)<br><sub>21-Subtype taxonomy; colors grouped by donor chemistry</sub>",
            'x': 0.5,
            'xanchor': 'center',
        },
        height=max(420, 60 + 26 * len(labels)),  # ~26 px per bar
        xaxis_title="Ligand count",
        yaxis_title=None,
        template="plotly_white",
        showlegend=False,
        margin=dict(l=10, r=80, t=80, b=40),
        bargap=0.25,
    )
    fig.update_yaxes(tickfont=dict(size=11))

    return fig


def plot_selectivity_trend_pie(ligands):
    """Creates a pie chart showing distribution of selectivity trends."""

    # Count selectivity trends
    trend_counts = {}
    for lig in ligands:
        trend_name = lig.selectivity_trend if lig.selectivity_trend else 'Unknown'
        trend_counts[trend_name] = trend_counts.get(trend_name, 0) + 1

    # Define custom order for trends (meaningful progression)
    trend_order = list(TREND_DISPLAY_ORDER)

    # Sort by custom order
    labels = []
    values = []
    for trend in trend_order:
        if trend in trend_counts:
            labels.append(trend)
            values.append(trend_counts[trend])

    # Add any trends not in predefined order
    for trend, count in trend_counts.items():
        if trend not in trend_order:
            labels.append(trend)
            values.append(count)

    # Define color scheme for trends
    color_map = {
        'Ascending': '#2E86AB',          # Blue (light→heavy preference)
        'Peaked': '#E76F51',             # Burnt orange (arch shape)
        'Flat': '#A23B72',               # Purple (no selectivity)
        'U_Shaped': '#457B9D',           # Steel blue (trough shape)
        'Reverse': '#F18F01',            # Orange (heavy→light preference)
        'Mixed': '#C73E1D',              # Red (conflicting trends)
        'Insufficient_Data': '#CCCCCC',  # Gray (not enough data)
        'Unknown': '#999999'             # Darker gray
    }
    colors = [color_map.get(label, '#999999') for label in labels]

    # Create pie chart
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.3,  # Donut chart
        hovertemplate='<b>%{label}</b><br>Count: %{value}<br>Percentage: %{percent}<extra></extra>',
        textposition='auto',
        textinfo='label+percent',
        marker=dict(
            colors=colors,
            line=dict(color='white', width=2)
        )
    )])

    fig.update_layout(
        title={
            'text': f"Selectivity Trend Distribution ({len(ligands)} total)",
            'x': 0.5,
            'xanchor': 'center'
        },
        height=400,
        template="plotly_white",
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="middle",
            y=0.5,
            xanchor="left",
            x=1.05
        )
    )

    return fig


def plot_lanthanide_distribution_donut(ligands):
    """Creates a donut chart showing distribution of measured lanthanides."""
    from collections import Counter

    # Count measurements for each lanthanide
    element_counts = Counter()
    for lig in ligands:
        # Count unique lanthanides across all conditions
        elements_in_ligand = set()
        for cond in lig.experimental_conditions:
            dist_ratios = cond.get('distribution_ratios', {})
            for elem in LANTHANIDES:
                if elem in dist_ratios:
                    elements_in_ligand.add(elem)

        # Add counts
        for elem in elements_in_ligand:
            element_counts[elem] += 1

    # Sort by atomic number
    labels = []
    values = []
    for elem in LANTHANIDES:
        if elem in element_counts:
            labels.append(elem)
            values.append(element_counts[elem])

    # Color scheme - gradient from light to heavy REE
    # Light REE (La-Nd): Blues
    # Middle REE (Pm-Tb): Greens
    # Heavy REE (Dy-Lu): Reds/Oranges
    color_map = {
        'La': '#E3F2FD', 'Ce': '#BBDEFB', 'Pr': '#90CAF9', 'Nd': '#64B5F6',
        'Pm': '#B9F6CA', 'Sm': '#69F0AE', 'Eu': '#00E676',
        'Gd': '#C8E6C9', 'Tb': '#81C784',
        'Dy': '#FFE082', 'Ho': '#FFD54F', 'Er': '#FFCA28',
        'Tm': '#FFB74D', 'Yb': '#FF9800', 'Lu': '#F57C00'
    }
    colors = [color_map.get(label, '#999999') for label in labels]

    # Create donut chart
    fig = go.Figure(data=[go.Pie(
        labels=labels,
        values=values,
        hole=0.4,  # Larger hole for donut
        hovertemplate='<b>%{label}</b><br>Ligands: %{value}<br>Percentage: %{percent}<extra></extra>',
        textposition='auto',
        textinfo='label+value',
        marker=dict(
            colors=colors,
            line=dict(color='white', width=2)
        )
    )])

    fig.update_layout(
        title={
            'text': f"Measured Lanthanides Distribution",
            'x': 0.5,
            'xanchor': 'center'
        },
        height=500,
        template="plotly_white",
        showlegend=False  # Too many elements for legend
    )

    return fig


def plot_element_distribution(df, element, x_axis_option):
    """
    Create interactive plot for single element distribution.

    Args:
        df: DataFrame from extract_element_data
        element: Element symbol
        x_axis_option: One of ["Ligand Types", "Solvents", "Solvent Categories", "Ligand Acidity",
                                "Acid Concentration", "Extractant Concentration", "Temperature"]

    Returns:
        Plotly figure
    """
    if df.empty:
        return _empty_plot_figure(f"No data available for {element}")

    # Load solvent categories if needed
    if x_axis_option == "Solvent Categories":
        solvent_to_category, category_descriptions = load_solvent_categories()

    # Determine x-axis data based on selection
    if x_axis_option == "Ligand Types":
        # Box plot grouped by Layer-1 Chemical_Family (13 broad categories).
        # The DataFrame also carries 'subtype' (Layer 2) for any future drill-down.
        plot_df = df.copy()
        x_column = 'chemical_family'
        x_label = "Ligand Family"
        plot_type = 'box'  # Box plot for categorical

    elif x_axis_option == "Solvents":
        # Keep raw data for box plot
        plot_df = df.copy()
        x_column = 'solvent'
        x_label = "Solvent"
        plot_type = 'box'  # Box plot for categorical

    elif x_axis_option == "Solvent Categories":
        # Keep raw data for box plot, add category column
        plot_df = df.copy()
        plot_df['category'] = plot_df['solvent'].apply(
            lambda s: get_solvent_category(s, solvent_to_category)
        )
        x_column = 'category'
        x_label = "Solvent Category"
        plot_type = 'box'  # Box plot for categorical

    elif x_axis_option == "Ligand Acidity":
        # Keep raw data for box plot
        plot_df = df.copy()
        x_column = 'acidity'
        x_label = "Ligand Acidity"
        plot_type = 'box'  # Box plot for categorical

    elif x_axis_option == "Acid Concentration":
        plot_df = df
        x_data = df['acid_conc']
        x_label = "Acid Concentration (M)"
        plot_type = 'scatter'

    elif x_axis_option == "Extractant Concentration":
        plot_df = df
        x_data = df['extractant_conc']
        x_label = "Extractant Concentration (M)"
        plot_type = 'scatter'

    elif x_axis_option == "Temperature":
        plot_df = df
        x_data = df['temperature']
        x_label = "Temperature (°C)"
        plot_type = 'scatter'

    # Statistics from values only; bounds are drawn on top as triangles
    plot_df, bound_df = split_bounds(plot_df)
    if plot_type == 'scatter':
        x_bounds = x_data.loc[bound_df.index]
        x_data = x_data.loc[plot_df.index]

    # Create figure
    fig = go.Figure()

    # Determine number of unique solvents for coloring (for scatter plots)
    unique_solvents = plot_df['solvent'].nunique() if 'solvent' in plot_df.columns else 1

    def get_category_color(category, x_axis_option):
        """
        Get color for a category based on x-axis option.

        Args:
            category: Category name (ligand type, solvent, solvent category, or acidity)
            x_axis_option: One of "Ligand Types", "Solvents", "Solvent Categories", "Ligand Acidity"

        Returns:
            Hex color string
        """
        if x_axis_option in ("Ligand Types", "Solvents",
                             "Solvent Categories", "Ligand Acidity"):
            return _category_color_map(x_axis_option).get(category, '#999999')
        else:
            return '#1f77b4'  # Default blue for unexpected cases

    if plot_type == 'box':
        # Box plot for categorical x-axis
        # Sort categories by median log_d of their values (descending); bound-only ones last
        categories_sorted = category_order(plot_df, bound_df, x_column, 'log_d')

        # Create one box trace per category
        for category in categories_sorted:
            mask = plot_df[x_column] == category
            category_data = plot_df[mask]
            n_bounds = int((bound_df[x_column] == category).sum())
            if category_data.empty:
                continue

            # Get color for this category
            category_color = get_category_color(category, x_axis_option)

            # Calculate statistics for hover
            n_measurements = len(category_data)
            median_val = category_data['log_d'].median()
            mean_val = category_data['log_d'].mean()
            q1_val = category_data['log_d'].quantile(0.25)
            q3_val = category_data['log_d'].quantile(0.75)
            min_val = category_data['log_d'].min()
            max_val = category_data['log_d'].max()

            # Get top 3 ligands for this category
            top_ligands = category_data.groupby('ligand_name')['log_d'].median().nlargest(3)
            ligand_examples = ', '.join(top_ligands.index.tolist())
            if len(ligand_examples) > 60:
                ligand_examples = ligand_examples[:60] + '...'

            # Category-specific info
            if x_axis_option == "Ligand Types":
                n_ligands = category_data['ligand_name'].nunique()
                extra_info = f"Unique ligands: {n_ligands}"
            elif x_axis_option == "Solvents":
                extra_info = f"Measurements: {n_measurements}"
            elif x_axis_option == "Solvent Categories":
                n_solvents = category_data['solvent'].nunique()
                extra_info = f"Unique solvents: {n_solvents}"
            elif x_axis_option == "Ligand Acidity":
                n_ligands = category_data['ligand_name'].nunique()
                extra_info = f"Unique ligands: {n_ligands}"

            # Custom hover template with statistics
            hover_template = (
                f"<b>{category}</b><br>"
                f"Median: {median_val:.2f}<br>"
                f"Mean: {mean_val:.2f}<br>"
                f"Q1: {q1_val:.2f}, Q3: {q3_val:.2f}<br>"
                f"Min: {min_val:.2f}, Max: {max_val:.2f}<br>"
                f"Count: {n_measurements}<br>"
                + (f"Bounds, not in these statistics: {n_bounds}<br>" if n_bounds else "")
                + f"{extra_info}<br>"
                f"Top ligands: {ligand_examples}"
                "<extra></extra>"
            )

            # Add violin plot FIRST (renders underneath)
            fig.add_trace(go.Violin(
                y=category_data['log_d'],
                name=category,
                line_color=category_color,
                fillcolor=category_color,
                opacity=0.3,  # Semi-transparent
                box_visible=False,  # Don't show violin's internal box
                meanline_visible=False,  # Don't show violin's mean line
                points=False,  # Don't show individual points (box plot will show outliers)
                scalemode='width',  # All violins have same width
                width=0.8,  # Slightly narrower than default
                showlegend=False,
                hoverinfo='skip'  # No hover on violin (box hover is better)
            ))

            # Add box plot on top with matching color
            fig.add_trace(go.Box(
                y=category_data['log_d'],
                name=category,
                marker=dict(
                    color=category_color,
                    line=dict(width=1.5, color='white')  # White border for visibility
                ),
                line=dict(color=category_color, width=2),  # Box outline
                boxmean=True,  # Show mean line
                boxpoints='outliers',  # Show only outliers
                hovertemplate=hover_template,
                showlegend=False,
                fillcolor='rgba(255,255,255,0.8)'  # Semi-transparent white fill for box
            ))

        # Add detailed hover for outlier points
        for category in categories_sorted:
            mask = plot_df[x_column] == category
            category_data = plot_df[mask]

            # Get color for this category
            category_color = get_category_color(category, x_axis_option)

            # Identify outliers using IQR method
            Q1 = category_data['log_d'].quantile(0.25)
            Q3 = category_data['log_d'].quantile(0.75)
            IQR = Q3 - Q1
            lower_fence = Q1 - 1.5 * IQR
            upper_fence = Q3 + 1.5 * IQR

            outliers = category_data[
                (category_data['log_d'] < lower_fence) |
                (category_data['log_d'] > upper_fence)
            ]

            if len(outliers) > 0:
                # Create detailed hover for each outlier
                hover_text_outliers = []
                for _, row in outliers.iterrows():
                    hover = (
                        f"<b>{row['ligand_name']}</b><br>"
                        f"log(D[{element}]): {row['log_d']:.2f}<br>"
                        f"D[{element}]: {row['d_value']:.2f}<br>"
                        f"Solvent: {row['solvent']}<br>"
                        f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"
                        f"Extractant Conc: {row['extractant_conc']:.3f} M<br>"
                        f"Temperature: {row['temperature']:.1f} °C"
                    )
                    hover_text_outliers.append(hover)

                # Add invisible scatter for outlier hover
                fig.add_trace(go.Scatter(
                    x=[category] * len(outliers),
                    y=outliers['log_d'],
                    mode='markers',
                    marker=dict(size=8, color=category_color, opacity=0),  # Invisible, matches category
                    hovertemplate='%{text}<extra></extra>',
                    text=hover_text_outliers,
                    showlegend=False,
                    hoverinfo='text'
                ))

        # Bounds on top, near-black so they stand out over the box
        for category in categories_sorted:
            cat_bounds = bound_df[bound_df[x_column] == category]
            if not cat_bounds.empty:
                fig.add_trace(bound_trace(
                    [category] * len(cat_bounds), cat_bounds['log_d'], cat_bounds['censor'],
                    [_element_bound_hover(r, element) for _, r in cat_bounds.iterrows()],
                    BOUND_COLOR))
    else:
        # Scatter plot for continuous x-axis
        # Generate hover text for scatter plots
        hover_text = [
            f"<b>{row['ligand_name']}</b><br>" +
            f"log(D[{element}]): {row['log_d']:.2f}<br>" +
            f"D[{element}]: {row['d_value']:.2f}<br>" +
            f"Solvent: {row['solvent']}<br>" +
            f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>" +
            f"Extractant Conc: {row['extractant_conc']:.3f} M<br>" +
            f"Temperature: {row['temperature']:.1f} °C"
            for _, row in plot_df.iterrows()
        ]

        # Color by solvent if multiple solvents present
        if unique_solvents <= 10:
            # Color by solvent (up to 10 unique solvents for readability)
            solvent_colors = {
                'n-dodecane': '#1f77b4',
                'toluene': '#ff7f0e',
                'chloroform': '#2ca02c',
                'kerosene': '#d62728',
                '1-octanol': '#9467bd',
                'Dichloromethane': '#8c564b',
                'meta-nitrobenzotrifluoride': '#e377c2',
                'benzene': '#7f7f7f',
                'nitrobenzene': '#bcbd22',
                'Unknown': '#17becf'
            }

            for solvent in plot_df['solvent'].unique():
                mask = plot_df['solvent'] == solvent
                hover_subset = [hover_text[i] for i, m in enumerate(mask) if m]

                fig.add_trace(go.Scatter(
                    x=x_data[mask],
                    y=plot_df.loc[mask, 'log_d'],
                    mode='markers',
                    name=solvent[:30],  # Truncate long solvent names
                    marker=dict(
                        size=10,
                        color=solvent_colors.get(solvent, '#7f7f7f'),
                        line=dict(width=1, color='white')
                    ),
                    hovertemplate='%{text}<extra></extra>',
                    text=hover_subset
                ))
        else:
            # Too many solvents - use single color
            fig.add_trace(go.Scatter(
                x=x_data,
                y=plot_df['log_d'],
                mode='markers',
                marker=dict(
                    size=10,
                    color='#1f77b4',
                    line=dict(width=1, color='black')
                ),
                hovertemplate='%{text}<extra></extra>',
                text=hover_text
            ))
        if not bound_df.empty:
            fig.add_trace(bound_trace(x_bounds, bound_df['log_d'], bound_df['censor'],
                                      [_element_bound_hover(r, element) for _, r in bound_df.iterrows()],
                                      BOUND_COLOR, showlegend=True))

    # Update layout
    # Show legend only for solvent-colored scatter plots
    show_legend = (plot_type == 'scatter' and unique_solvents <= 10)

    font_sizes = get_font_sizes()
    fig.update_layout(
        title=dict(
            text=f"Distribution Coefficient for {element}",
            font=dict(size=font_sizes['title'], family="Arial, sans-serif", color="#2C3E50")
        ),
        xaxis=dict(
            title=dict(
                text=x_label,
                font=dict(size=font_sizes['axis_title'], family="Arial, sans-serif")
            ),
            tickfont=dict(size=font_sizes['tick']),
            gridcolor='rgba(200,200,200,0.3)'
        ),
        yaxis=dict(
            title=dict(
                text=f"log(D[{element}])",
                font=dict(size=font_sizes['axis_title'], family="Arial, sans-serif")
            ),
            tickfont=dict(size=font_sizes['tick']),
            gridcolor='rgba(200,200,200,0.3)'
        ),
        template="plotly_white",
        height=500,
        hovermode='closest',
        plot_bgcolor='rgba(250,250,250,0.5)',
        hoverlabel=dict(
            font_size=font_sizes['hover'],
            font_family="Arial"
        ),
        showlegend=show_legend,
        legend=dict(
            title=dict(text="Solvent", font=dict(size=font_sizes['legend'])),
            font=dict(size=font_sizes['legend']),
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01
        ) if show_legend else None
    )

    # Rotate x-axis labels if categorical (ligand types, solvents, solvent categories, or ligand acidity)
    if x_axis_option in ["Ligand Types", "Solvents", "Solvent Categories", "Ligand Acidity"]:
        fig.update_xaxes(tickangle=-45)

        # For box plots, enforce category order (sorted by median descending)
        if plot_type == 'box':
            fig.update_xaxes(
                categoryorder='array',
                categoryarray=categories_sorted
            )

    return fig


def plot_separation_factor_distribution(df, element1, element2, x_axis_option):
    """
    Create interactive plot for separation factor distribution.

    Args:
        df: DataFrame from extract_separation_factor_data
        element1: First element symbol (e.g., 'Eu')
        element2: Second element symbol (e.g., 'La')
        x_axis_option: One of 7 x-axis options

    Returns:
        Plotly figure with logSF distribution
    """
    # Handle empty dataframe
    if df.empty:
        return _empty_plot_figure(
            f"No paired data available for {element1} and {element2}",
            title=f"Separation Factor: {element1} vs {element2}",
            hide_axes=True,
        )

    # Determine x-axis column and plot type based on option
    x_column = None
    plot_type = None  # 'box' for categorical, 'scatter' for continuous

    if x_axis_option == "Ligand Types":
        x_column = 'chemical_family'
        plot_type = 'box'
    elif x_axis_option == "Solvents":
        x_column = 'solvent'
        plot_type = 'box'
    elif x_axis_option == "Solvent Categories":
        # Add solvent category column
        df = df.copy()
        solvent_to_category, _ = load_solvent_categories()
        df['solvent_category'] = df['solvent'].apply(
            lambda s: get_solvent_category(s, solvent_to_category)
        )
        x_column = 'solvent_category'
        plot_type = 'box'
    elif x_axis_option == "Ligand Acidity":
        x_column = 'acidity'
        plot_type = 'box'
    elif x_axis_option == "Acid Concentration":
        x_column = 'acid_conc'
        plot_type = 'scatter'
    elif x_axis_option == "Extractant Concentration":
        x_column = 'extractant_conc'
        plot_type = 'scatter'
    elif x_axis_option == "Temperature":
        x_column = 'temperature'
        plot_type = 'scatter'

    # Filter out rows with missing x-values
    plot_df = df[df[x_column].notna()].copy()

    if plot_df.empty:
        return _empty_plot_figure(
            f"No valid {x_axis_option.lower()} data available",
            title=f"Separation Factor: {element1} vs {element2}",
        )

    # Statistics from values only; bounds are drawn on top as triangles
    plot_df, bound_df = split_bounds(plot_df)

    # Create figure
    fig = go.Figure()

    # Box + Violin plots for categorical x-axis
    if plot_type == 'box':
        # Get unique categories
        categories = plot_df[x_column].unique()

        # Sort categories by median logSF of their values (descending); bound-only ones last
        categories = category_order(plot_df, bound_df, x_column, 'log_sf')

        # Determine color scheme
        color_map = _category_color_map(x_axis_option)

        # Get font sizes for responsive plots
        font_sizes = get_font_sizes()

        # Add violin plots (background layer)
        for category in categories:
            category_data = plot_df[plot_df[x_column] == category]
            color = color_map.get(category, '#CCCCCC')
            n_bounds = int((bound_df[x_column] == category).sum())
            if category_data.empty:
                continue

            # Build hover text for violin (shows stats + top ligands)
            log_sf_values = category_data['log_sf'].values
            median_val = np.median(log_sf_values)
            mean_val = np.mean(log_sf_values)
            q1 = np.percentile(log_sf_values, 25)
            q3 = np.percentile(log_sf_values, 75)
            min_val = np.min(log_sf_values)
            max_val = np.max(log_sf_values)

            # Get top ligands by logSF
            top_ligands = category_data.nlargest(3, 'log_sf')['ligand_name'].tolist()
            ligands_str = '<br>'.join([f"  • {lig}" for lig in top_ligands[:3]])

            hover_text = (
                f"<b>{category}</b><br>"
                f"N = {len(category_data)}<br>"
                + (f"Bounds, not in these statistics: {n_bounds}<br>" if n_bounds else "")
                + f"<br>"
                f"Median logSF: {median_val:.2f}<br>"
                f"Mean logSF: {mean_val:.2f}<br>"
                f"Q1: {q1:.2f}<br>"
                f"Q3: {q3:.2f}<br>"
                f"Min: {min_val:.2f}<br>"
                f"Max: {max_val:.2f}<br>"
                f"<br>"
                f"Top ligands:<br>{ligands_str}"
            )

            fig.add_trace(go.Violin(
                x=[category] * len(category_data),
                y=category_data['log_sf'],
                name=category,
                fillcolor=color,
                opacity=0.4,
                line_color=color,
                line_width=2,
                showlegend=False,
                scalemode='width',
                width=0.8,
                hoverinfo='text',
                hovertext=[hover_text] * len(category_data),
                hoverlabel=dict(font=dict(size=font_sizes['hover']))
            ))

        # Add box plots (foreground layer)
        for category in categories:
            category_data = plot_df[plot_df[x_column] == category]
            color = color_map.get(category, '#CCCCCC')
            if category_data.empty:
                continue

            # Build detailed hover text for individual points
            hover_texts = []
            for _, row in category_data.iterrows():
                hover_text = (
                    f"<b>{row['ligand_name']}</b><br>"
                    f"logSF[{element1}/{element2}]: {row['log_sf']:.2f}<br>"
                    f"logD[{element1}]: {row['log_d_elem1']:.2f}<br>"
                    f"logD[{element2}]: {row['log_d_elem2']:.2f}<br>"
                    f"<br>"
                    f"Solvent: {row['solvent']}<br>"
                    f"Acid: {row['acid_name']}"
                )
                if row['acid_conc'] is not None:
                    hover_text += f" ({row['acid_conc']:.2f} M)"
                if row['extractant_conc'] is not None:
                    hover_text += f"<br>Extractant conc: {row['extractant_conc']:.2f} M"
                if row['temperature'] is not None:
                    hover_text += f"<br>Temperature: {row['temperature']:.1f} °C"

                hover_texts.append(hover_text)

            fig.add_trace(go.Box(
                x=[category] * len(category_data),
                y=category_data['log_sf'],
                name=category,
                marker_color=color,
                fillcolor='rgba(255,255,255,0.8)',
                line=dict(color=color, width=3),
                boxmean='sd',
                showlegend=False,
                hoverinfo='text',
                hovertext=hover_texts,
                hoverlabel=dict(font=dict(size=font_sizes['hover'])),
                marker=dict(
                    size=12,
                    line=dict(color=color, width=2)
                )
            ))

        # Bounds on top, near-black so they stand out over the box
        for category in categories:
            cat_bounds = bound_df[bound_df[x_column] == category]
            if not cat_bounds.empty:
                fig.add_trace(bound_trace(
                    [category] * len(cat_bounds), cat_bounds['log_sf'], cat_bounds['censor'],
                    [_sf_bound_hover(r, element1, element2) for _, r in cat_bounds.iterrows()],
                    BOUND_COLOR))

    # Scatter plots for continuous x-axis
    elif plot_type == 'scatter':
        # Determine if we should color by solvent
        unique_solvents = plot_df['solvent'].nunique()

        # Get font sizes for responsive plots
        font_sizes = get_font_sizes()

        if unique_solvents <= 10:
            # Color by solvent
            for solvent in plot_df['solvent'].unique():
                solvent_data = plot_df[plot_df['solvent'] == solvent]
                color = SOLVENT_COLORS.get(solvent, '#CCCCCC')

                # Build hover text
                hover_texts = []
                for _, row in solvent_data.iterrows():
                    hover_text = (
                        f"<b>{row['ligand_name']}</b><br>"
                        f"logSF[{element1}/{element2}]: {row['log_sf']:.2f}<br>"
                        f"logD[{element1}]: {row['log_d_elem1']:.2f}<br>"
                        f"logD[{element2}]: {row['log_d_elem2']:.2f}<br>"
                        f"<br>"
                        f"Solvent: {row['solvent']}<br>"
                    )
                    if x_axis_option == "Acid Concentration":
                        hover_text += f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"
                        if row['extractant_conc'] is not None:
                            hover_text += f"Extractant conc: {row['extractant_conc']:.2f} M<br>"
                    elif x_axis_option == "Extractant Concentration":
                        hover_text += f"Extractant conc: {row['extractant_conc']:.2f} M<br>"
                        if row['acid_conc'] is not None:
                            hover_text += f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"
                    elif x_axis_option == "Temperature":
                        hover_text += f"Temperature: {row['temperature']:.1f} °C<br>"
                        if row['acid_conc'] is not None:
                            hover_text += f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"

                    hover_texts.append(hover_text)

                fig.add_trace(go.Scatter(
                    x=solvent_data[x_column],
                    y=solvent_data['log_sf'],
                    mode='markers',
                    name=solvent,
                    marker=dict(
                        size=16,
                        color=color,
                        line=dict(color='white', width=2),
                        opacity=0.8
                    ),
                    hoverinfo='text',
                    hovertext=hover_texts,
                    hoverlabel=dict(font=dict(size=font_sizes['hover']))
                ))
        else:
            # Too many solvents, use single color
            hover_texts = []
            for _, row in plot_df.iterrows():
                hover_text = (
                    f"<b>{row['ligand_name']}</b><br>"
                    f"logSF[{element1}/{element2}]: {row['log_sf']:.2f}<br>"
                    f"logD[{element1}]: {row['log_d_elem1']:.2f}<br>"
                    f"logD[{element2}]: {row['log_d_elem2']:.2f}<br>"
                    f"<br>"
                    f"Solvent: {row['solvent']}<br>"
                )
                if x_axis_option == "Acid Concentration":
                    hover_text += f"Acid: {row['acid_name']} ({row['acid_conc']:.2f} M)<br>"
                elif x_axis_option == "Extractant Concentration":
                    hover_text += f"Extractant conc: {row['extractant_conc']:.2f} M<br>"
                elif x_axis_option == "Temperature":
                    hover_text += f"Temperature: {row['temperature']:.1f} °C<br>"

                hover_texts.append(hover_text)

            fig.add_trace(go.Scatter(
                x=plot_df[x_column],
                y=plot_df['log_sf'],
                mode='markers',
                name=element1,
                marker=dict(
                    size=16,
                    color='#1f77b4',
                    line=dict(color='white', width=2),
                    opacity=0.8
                ),
                hoverinfo='text',
                hovertext=hover_texts,
                hoverlabel=dict(font=dict(size=font_sizes['hover'])),
                showlegend=False
            ))
        if not bound_df.empty:
            fig.add_trace(bound_trace(bound_df[x_column], bound_df['log_sf'], bound_df['censor'],
                                      [_sf_bound_hover(r, element1, element2) for _, r in bound_df.iterrows()],
                                      BOUND_COLOR, showlegend=True))

    # Add reference line at y=0 (no selectivity)
    # Note: font_sizes is already available from the if/elif blocks above
    if 'font_sizes' not in locals():
        font_sizes = get_font_sizes()

    fig.add_hline(
        y=0,
        line_dash="dash",
        line_color="gray",
        line_width=2,
        annotation_text="No selectivity (logSF=0)",
        annotation_position="right",
        annotation_font_size=font_sizes['annotation']
    )

    # Update layout
    x_axis_labels = {
        "Ligand Types": "Extractant Type",
        "Solvents": "Organic Solvent",
        "Solvent Categories": "Solvent Category",
        "Ligand Acidity": "Ligand Acidity (Extraction Mechanism)",
        "Acid Concentration": "Acid Concentration (M)",
        "Extractant Concentration": "Extractant Concentration (M)",
        "Temperature": "Temperature (°C)"
    }

    # Ensure font_sizes is available (should be set from if/elif blocks above)
    if 'font_sizes' not in locals():
        font_sizes = get_font_sizes()

    fig.update_layout(
        title=dict(
            text=f"Separation Factor: {element1} vs {element2}",
            font=dict(size=font_sizes['title'], family="Arial, sans-serif", color="#2C3E50"),
            x=0.5,
            xanchor='center'
        ),
        xaxis_title=dict(
            text=x_axis_labels[x_axis_option],
            font=dict(size=font_sizes['axis_title'], family="Arial, sans-serif")
        ),
        yaxis_title=dict(
            text=f"logSF[{element1}/{element2}] = logD[{element1}] - logD[{element2}]",
            font=dict(size=font_sizes['axis_title'], family="Arial, sans-serif")
        ),
        xaxis=dict(
            tickfont=dict(size=font_sizes['tick']),
            gridcolor='rgba(200,200,200,0.3)'
        ),
        yaxis=dict(
            tickfont=dict(size=font_sizes['tick']),
            gridcolor='rgba(200,200,200,0.3)'
        ),
        height=500,
        template="plotly_white",
        hovermode='closest',
        plot_bgcolor='rgba(250,250,250,0.5)',
        showlegend=(plot_type == 'scatter' and unique_solvents <= 10)
    )

    # Rotate x-axis labels if categorical
    if plot_type == 'box':
        fig.update_xaxes(tickangle=-45)

        # Enforce category order (sorted by median descending)
        fig.update_xaxes(
            categoryorder='array',
            categoryarray=categories
        )

    return fig


def plot_eu_la_scatter(ligands):
    """Creates scatter plots for Eu/La separation factors."""
    from plotly.subplots import make_subplots

    # Extract Eu and La data
    data_points = []
    for lig in ligands:
        # Check if ligand has both Eu and La
        tested = lig.tested_elements
        if 'Eu' not in tested or 'La' not in tested:
            continue

        # Get aggregated D-values
        try:
            d_values = lig.get_aggregated_d_values(aggregation='median', lanthanides_only=True)

            if 'Eu' in d_values and 'La' in d_values:
                log_d_eu = d_values['Eu']['log_value']
                log_d_la = d_values['La']['log_value']

                if log_d_eu is not None and log_d_la is not None:
                    # Calculate separation factor: log(SF[La/Eu]) = log(D[La]) - log(D[Eu])
                    # Positive = La-selective, Negative = Eu-selective
                    log_sf_la_eu = log_d_la - log_d_eu

                    data_points.append({
                        'name': lig.name,
                        'smiles': lig.smiles if lig.smiles else 'N/A',
                        'log_d_eu': log_d_eu,
                        'log_d_la': log_d_la,
                        'log_sf': log_sf_la_eu,
                        'mechanism': lig.mechanism,
                        'structural_class': lig.structural_class
                    })
        except Exception as e:
            continue

    if not data_points:
        # Return empty figure if no data
        fig = go.Figure()
        fig.add_annotation(
            text="No Eu/La data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig

    df = pd.DataFrame(data_points)

    # Calculate common axis ranges for consistency
    x_min = min(df['log_d_eu'].min(), df['log_d_la'].min())
    x_max = max(df['log_d_eu'].max(), df['log_d_la'].max())
    x_margin = (x_max - x_min) * 0.05
    x_range = [x_min - x_margin, x_max + x_margin]

    y_min = df['log_sf'].min()
    y_max = df['log_sf'].max()
    y_margin = (y_max - y_min) * 0.05
    y_range = [y_min - y_margin, y_max + y_margin]

    # Create subplots
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=('log(D[Eu]) vs log(SF[La/Eu])',
                       'log(D[La]) vs log(SF[La/Eu])'),
        horizontal_spacing=0.12
    )

    # Plot 1: Eu
    hover_text_eu = [
        f"<b>{row['name']}</b><br>" +
        f"SMILES: {row['smiles'][:50]}{'...' if len(row['smiles']) > 50 else ''}<br>" +
        f"log(D[Eu]): {row['log_d_eu']:.2f}<br>" +
        f"log(D[La]): {row['log_d_la']:.2f}<br>" +
        f"log(SF[La/Eu]): {row['log_sf']:.2f}<br>" +
        f"SF[La/Eu]: {10**row['log_sf']:.2f}<br>" +
        f"Mechanism: {row['mechanism']}"
        for _, row in df.iterrows()
    ]

    fig.add_trace(
        go.Scatter(
            x=df['log_d_eu'],
            y=df['log_sf'],
            mode='markers',
            marker=dict(size=8, color='#1f77b4', line=dict(width=1, color='black')),
            hovertemplate='%{text}<extra></extra>',
            text=hover_text_eu,
            showlegend=False
        ),
        row=1, col=1
    )

    # Plot 2: La
    hover_text_la = [
        f"<b>{row['name']}</b><br>" +
        f"SMILES: {row['smiles'][:50]}{'...' if len(row['smiles']) > 50 else ''}<br>" +
        f"log(D[Eu]): {row['log_d_eu']:.2f}<br>" +
        f"log(D[La]): {row['log_d_la']:.2f}<br>" +
        f"log(SF[La/Eu]): {row['log_sf']:.2f}<br>" +
        f"SF[La/Eu]: {10**row['log_sf']:.2f}<br>" +
        f"Mechanism: {row['mechanism']}"
        for _, row in df.iterrows()
    ]

    fig.add_trace(
        go.Scatter(
            x=df['log_d_la'],
            y=df['log_sf'],
            mode='markers',
            marker=dict(size=8, color='#1f77b4', line=dict(width=1, color='black')),
            hovertemplate='%{text}<extra></extra>',
            text=hover_text_la,
            showlegend=False
        ),
        row=1, col=2
    )

    # Add reference lines at y=0
    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5,
                  annotation_text="No selectivity", row=1, col=1)
    fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5,
                  annotation_text="No selectivity", row=1, col=2)

    # Update axes with common ranges
    fig.update_xaxes(title_text="log(D[Eu])", row=1, col=1, showgrid=True, range=x_range)
    fig.update_xaxes(title_text="log(D[La])", row=1, col=2, showgrid=True, range=x_range)
    fig.update_yaxes(title_text="log(SF[La/Eu])", row=1, col=1, showgrid=True, range=y_range)
    fig.update_yaxes(title_text="log(SF[La/Eu])", row=1, col=2, showgrid=True, range=y_range)

    # Update layout
    fig.update_layout(
        template="plotly_white",
        height=400,
        hovermode='closest'
    )

    return fig


#=============================================================================
# AI Assistant Tab
#=============================================================================

# Suggested questions shown as quick-start buttons
_EXAMPLE_QUESTIONS = [
    "What are the best ligands for separating Eu from La?",
    "Show me DGA ligands with ascending selectivity trend",
    "Which ligands have data for all 15 lanthanides?",
    "Compare TODGA and TEHDGA",
    "What is the database statistics summary?",
    "Find organophosphorus ligands with reverse selectivity",
    "Which ligands have the highest logD for Nd?",
    "What ligands are available from the OECD-NEA IDEAL database?",
]

# Curated OpenRouter free-tier models offered in free mode. Any other free
# model ID can be typed in, but these are the ones checked to work with the
# assistant's context-injection prompt. "openrouter/free" leads because it is a
# router rather than one model — it picks whichever free models are currently
# up, so it cannot go stale the way a pinned ID does.
_FREE_MODELS = [
    "openrouter/free",
    "nvidia/nemotron-3-ultra-550b-a55b:free",
    "nvidia/nemotron-3-nano-omni-30b-a3b-reasoning:free",
    "minimax/minimax-m3.1-flash:free",
]
_CUSTOM_FREE_LABEL = "Other :free model…"

_MODE_LABELS = {
    "free": "🎁 Free models",
    "user": "🔑 Your own API key",
    "ollama": "🖥️ Local (Ollama)",
}


def _free_key() -> str:
    """The app owner's shared OpenRouter key, if one is configured."""
    return get_config("OPENROUTER_FREE_API_KEY")


def _ollama_offered() -> bool:
    """
    Whether to show the local-Ollama mode.

    Probed once per session: on Streamlit Cloud nothing listens on localhost,
    so the connection is refused immediately and cloud visitors never see a
    dead option. Set RELQUEST_ENABLE_OLLAMA=1 to show it regardless.
    """
    if get_config("RELQUEST_ENABLE_OLLAMA"):
        return True
    if "ollama_reachable" not in st.session_state:
        st.session_state.ollama_reachable = check_ollama_connection()
    return st.session_state.ollama_reachable


def _available_modes() -> list:
    """Backend modes offered in this deployment, best default first."""
    modes = []
    if _free_key():
        modes.append("free")
    modes.append("user")
    if _ollama_offered():
        modes.append("ollama")
    return modes


def _init_llm_settings():
    """Seed the LLM settings in session state on first render."""
    modes = _available_modes()
    if st.session_state.get("llm_mode") not in modes:
        st.session_state.llm_mode = modes[0]
    st.session_state.setdefault("llm_free_model", _FREE_MODELS[0])
    st.session_state.setdefault("llm_user_key", get_config("RELQUEST_API_KEY"))
    st.session_state.setdefault(
        "llm_user_base_url", get_config("RELQUEST_CLOUD_BASE_URL", OPENROUTER_BASE_URL))
    st.session_state.setdefault(
        "llm_user_model", get_config("RELQUEST_CLOUD_MODEL", "openai/gpt-4o-mini"))
    st.session_state.setdefault("llm_ollama_model", DEFAULT_MODEL)


@st.dialog("⚙️ AI Settings")
def _show_llm_settings_dialog():
    """Popup for choosing the LLM backend, model and (optionally) an API key."""
    modes = _available_modes()

    mode = st.radio(
        "How would you like to run the assistant?",
        options=modes,
        format_func=lambda m: _MODE_LABELS[m],
        index=modes.index(st.session_state.llm_mode),
        key="settings_mode_radio",
    )

    free_model = st.session_state.llm_free_model
    user_key = st.session_state.llm_user_key
    user_base_url = st.session_state.llm_user_base_url
    user_model = st.session_state.llm_user_model
    ollama_model = st.session_state.llm_ollama_model
    error = None

    if mode == "free":
        st.markdown("#### Free models")
        st.info(
            "These run on the app owner's OpenRouter key — no signup needed. "
            "The free tier has a small daily quota shared by everyone using "
            "this app, so switch to your own key if you hit a rate limit.",
            icon="🎁",
        )
        options = _FREE_MODELS + [_CUSTOM_FREE_LABEL]
        preset = free_model if free_model in _FREE_MODELS else _CUSTOM_FREE_LABEL
        choice = st.selectbox(
            "Model", options, index=options.index(preset), key="settings_free_model")
        if choice == _CUSTOM_FREE_LABEL:
            custom = st.text_input(
                "OpenRouter model ID",
                value="" if free_model in _FREE_MODELS else free_model,
                placeholder="e.g. qwen/qwen-2.5-72b-instruct:free",
                help="Must end in `:free` — the shared key is limited to free models.",
                key="settings_free_model_custom",
            )
            st.markdown(
                "🔗 [Browse free models on OpenRouter](https://openrouter.ai/models?q=free)")
            if custom and not is_free_model(custom):
                error = (f"`{custom}` is not a free model. Free-tier model IDs end "
                         "in `:free`. To use a paid model, choose **Your own API key**.")
            free_model = custom.strip() if custom else free_model
        else:
            free_model = choice

    elif mode == "user":
        st.markdown("#### Your own API key")
        st.info(
            "Works with any OpenAI-compatible provider — OpenRouter, OpenAI, "
            "DeepSeek, Groq, or a local proxy. Your key is used only for this "
            "session and is never stored on the server.",
            icon="🔑",
        )
        user_key = st.text_input(
            "API key", type="password", value=user_key,
            help="Free OpenRouter keys: https://openrouter.ai/keys",
            key="settings_user_key")
        user_base_url = st.text_input(
            "Base URL", value=user_base_url, key="settings_user_base_url")
        user_model = st.text_input(
            "Model", value=user_model, key="settings_user_model")

        with st.expander("📋 Provider examples"):
            st.markdown(
                "| Provider | Base URL | Example model |\n"
                "|---|---|---|\n"
                "| Anthropic (Claude) | `https://api.anthropic.com/v1/` | `claude-opus-5` |\n"
                "| OpenRouter | `https://openrouter.ai/api/v1` | `anthropic/claude-opus-5` |\n"
                "| OpenAI | `https://api.openai.com/v1` | `gpt-4o-mini` |\n"
                "| DeepSeek | `https://api.deepseek.com/v1` | `deepseek-chat` |\n"
                "| Groq | `https://api.groq.com/openai/v1` | `llama-3.3-70b-versatile` |\n"
                "\n"
                "**Claude** — an Anthropic key works directly against the base URL "
                "above (keep the trailing slash), or via OpenRouter with an "
                "`anthropic/` prefix. `claude-sonnet-5` is the cheaper option; "
                "`claude-haiku-4-5` cheaper still."
            )
        if not user_key:
            st.warning("Enter an API key to enable the assistant.", icon="⚠️")

    else:
        st.markdown("#### Local Ollama")
        if check_ollama_connection():
            st.success("Ollama connected", icon="✅")
            models = list_available_models()
        else:
            st.error("Ollama is not reachable at this address.", icon="❌")
            st.markdown("```bash\nollama serve\n# or\ndocker compose up -d ollama\n```")
            models = []
        options = models or [DEFAULT_MODEL]
        if ollama_model not in options:
            options = [ollama_model] + options
        ollama_model = st.selectbox(
            "Model", options, index=options.index(ollama_model),
            help="Pull models with `ollama pull qwen2.5:7b` — qwen2.5:7b supports tool calling.",
            key="settings_ollama_model")

    if error:
        st.error(error, icon="🚫")

    if st.button("Save", type="primary", use_container_width=True, disabled=bool(error)):
        st.session_state.llm_mode = mode
        st.session_state.llm_free_model = free_model
        st.session_state.llm_user_key = user_key
        st.session_state.llm_user_base_url = user_base_url
        st.session_state.llm_user_model = user_model
        st.session_state.llm_ollama_model = ollama_model
        st.rerun()


def _configured_model() -> str:
    """The model the saved settings point at, for the status line."""
    return st.session_state[{
        "free": "llm_free_model",
        "user": "llm_user_model",
        "ollama": "llm_ollama_model",
    }[st.session_state.llm_mode]]


def _build_llm_backend():
    """
    Construct the backend for the saved settings.

    Returns (backend, model, base_url, api_key). backend is None when the
    configuration is incomplete or unusable; the caller reports why.
    """
    mode = st.session_state.llm_mode

    if mode == "free":
        key, model = _free_key(), st.session_state.llm_free_model
        if not is_free_model(model):
            st.error(
                f"`{model}` is not a free-tier model, so it cannot use the shared key. "
                "Open ⚙️ Settings to pick a `:free` model or supply your own key.",
                icon="🚫")
            return None, model, OPENROUTER_BASE_URL, ""
        return _make_cloud_backend(key, model, OPENROUTER_BASE_URL)

    if mode == "user":
        key = st.session_state.llm_user_key
        model = st.session_state.llm_user_model
        base_url = st.session_state.llm_user_base_url
        if not key:
            st.warning(
                "Open ⚙️ Settings and enter an API key to start chatting.", icon="🔑")
            return None, model, base_url, ""
        return _make_cloud_backend(key, model, base_url)

    model = st.session_state.llm_ollama_model
    if not check_ollama_connection():
        st.error(
            "Ollama is not reachable. Start it with `ollama serve`, or pick another "
            "backend in ⚙️ Settings.", icon="❌")
        return None, model, OLLAMA_BASE_URL, ""
    return OllamaBackend(OLLAMA_BASE_URL, model), model, OLLAMA_BASE_URL, ""


def _make_cloud_backend(api_key: str, model: str, base_url: str):
    """Instantiate the OpenAI-compatible backend, reporting a missing SDK."""
    try:
        return OpenAICompatibleBackend(api_key, model, base_url), model, base_url, api_key
    except ImportError:
        st.error(
            "The `openai` package is required for cloud backends. "
            "Install it with: `pip install openai>=1.0.0`", icon="📦")
        return None, model, base_url, api_key


def _render_ai_assistant_tab(db: LigandDatabase):
    """Render the AI Assistant chat interface."""

    st.subheader("🤖 REEL-QUEST AI Assistant")
    st.markdown(
        "Ask questions about the ligand database in natural language. "
        "Answers are grounded in real data — not hallucinated."
    )

    if not AGENT_AVAILABLE:
        st.error(
            "The `requests` package is required for the AI assistant. "
            "Install it with: `pip install requests`"
        )
        return

    # ── Backend status + settings ────────────────────────────────────────────
    _init_llm_settings()
    mode = st.session_state.llm_mode

    col_status, col_button = st.columns([5, 1])
    with col_status:
        st.caption(f"{_MODE_LABELS[mode]} · `{_configured_model()}`")
    with col_button:
        if st.button("⚙️ Settings", use_container_width=True):
            _show_llm_settings_dialog()

    # Built after the status row so any problem is reported beneath it.
    backend, selected_model, active_base_url, active_api_key = _build_llm_backend()
    is_connected = backend is not None

    if is_connected and mode == "free":
        st.caption(
            "Free models are shared by everyone using this app and have a small "
            "daily quota. Hitting a rate limit? Add your own key in ⚙️ Settings — "
            "[free OpenRouter keys](https://openrouter.ai/keys)."
        )

    st.markdown("---")

    # ── Session state ────────────────────────────────────────────────────────
    if "agent_messages" not in st.session_state:
        st.session_state.agent_messages = []
    if "agent_instance" not in st.session_state:
        st.session_state.agent_instance = None
    if "agent_backend_key" not in st.session_state:
        st.session_state.agent_backend_key = None

    # Key that identifies the current backend configuration. It reads the
    # settings actually in force — an earlier version keyed on the environment,
    # so switching model or pasting a different key left the cached agent
    # untouched. The key itself is fingerprinted rather than parked verbatim in
    # session state.
    key_fp = (hashlib.sha256(active_api_key.encode()).hexdigest()[:12]
              if active_api_key else "")
    backend_key = (mode, selected_model, active_base_url, key_fp)

    # Re-create agent if backend configuration changed
    if st.session_state.agent_backend_key != backend_key:
        st.session_state.agent_instance = None
        st.session_state.agent_backend_key = backend_key

    if st.session_state.agent_instance is None and backend is not None:
        st.session_state.agent_instance = RELQuestAgentStreaming(
            db, backend=backend, model=selected_model
        )

    agent: RELQuestAgentStreaming = st.session_state.agent_instance

    # ── Example question buttons ─────────────────────────────────────────────
    if not st.session_state.agent_messages:
        st.markdown("**Quick start — click a question to ask it:**")
        cols = st.columns(2)
        for i, q in enumerate(_EXAMPLE_QUESTIONS):
            if cols[i % 2].button(q, key=f"eg_{i}", use_container_width=True):
                st.session_state._pending_question = q
                st.rerun()

    # ── Chat history display ─────────────────────────────────────────────────
    for msg in st.session_state.agent_messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])
            # Re-render any stored figures
            if msg.get("figures"):
                for fig in msg["figures"]:
                    st.plotly_chart(fig, use_container_width=True)

    # ── Handle pending question from example buttons ─────────────────────────
    if "agent_messages" in st.session_state and hasattr(st.session_state, "_pending_question"):
        pending = st.session_state._pending_question
        del st.session_state._pending_question
        if agent is not None:
            _handle_user_message(agent, pending)

    # ── Chat input ───────────────────────────────────────────────────────────
    user_input = st.chat_input(
        "Ask about ligands, trends, separation factors, mechanisms…",
        disabled=not is_connected,
    )

    if user_input and agent is not None:
        _handle_user_message(agent, user_input)

    # ── Clear conversation button ────────────────────────────────────────────
    if st.session_state.agent_messages and agent is not None:
        if st.button("🗑️ Clear conversation", use_container_width=False):
            st.session_state.agent_messages = []
            agent.reset()
            st.rerun()

    # ── Ollama model pull helper (only for local backend) ────────────────────
    if mode == "ollama" and is_connected:
        available_models = list_available_models()
        if not available_models:
            st.info(
                f"No models found. Pull the default model:\n"
                f"```bash\nollama pull {DEFAULT_MODEL}\n```"
            )
        elif selected_model not in available_models:
            with st.expander("⬇️ Pull this model"):
                st.markdown(f"Model `{selected_model}` is not yet downloaded.")
                if st.button(f"Pull {selected_model}", type="primary"):
                    from rel_quest.agent import pull_model
                    with st.spinner(f"Pulling {selected_model}… (this may take several minutes)"):
                        ok = pull_model(selected_model)
                    if ok:
                        st.success("Model pulled successfully. Refresh the page.")
                    else:
                        st.error("Pull failed. Run manually: `ollama pull " + selected_model + "`")


def _handle_user_message(agent: "RELQuestAgentStreaming", user_input: str):
    """Add user message to history and stream the assistant response."""
    st.session_state.agent_messages.append({"role": "user", "content": user_input})

    with st.chat_message("user"):
        st.markdown(user_input)

    with st.chat_message("assistant"):
        response_placeholder = st.empty()
        full_response = ""
        with st.spinner("Thinking…"):
            for chunk in agent.chat(user_input):
                full_response += chunk
                response_placeholder.markdown(full_response + "▌")
        response_placeholder.markdown(full_response)

        # Render any figures the agent generated
        if hasattr(agent, 'last_figures') and agent.last_figures:
            for fig in agent.last_figures:
                st.plotly_chart(fig, use_container_width=True)

    # Store message (and figure count for re-rendering)
    msg_entry = {"role": "assistant", "content": full_response}
    if hasattr(agent, 'last_figures') and agent.last_figures:
        msg_entry["figures"] = agent.last_figures
    st.session_state.agent_messages.append(msg_entry)


#=============================================================================
# Main Application
#=============================================================================

def _render_sidebar(db):
    """
    Render the sidebar: page navigation and every dashboard filter.

    Returns:
        tuple: (page, show_liu_only, search_text, extractant_types_filter,
        mechanisms, trends, sources, required_elements, logd_min, logd_max)
    """
    with st.sidebar:
        # ── Page navigation ──────────────────────────────────────────────────
        page = st.radio(
            "Navigation",
            ["📊 Database", "🤖 AI Assistant"],
            label_visibility="collapsed",
            horizontal=True,
        )
        st.markdown("---")

        # Show filters only on the Database page
        if page == "🤖 AI Assistant":
            st.info("Ask questions about the REEL-QUEST ligand database in natural language.")
        else:
            st.markdown('<div id="filters"></div>', unsafe_allow_html=True)
            st.header("Filters")

        # Search
        search_text = st.text_input("🔍 Search ligand name", "", help="Search by ligand name (case-insensitive)")

        # Type filter — uses Layer-1 Chemical_Family (13 broad categories: Diglycolamide,
        # Phenanthroline_Amide, Amide_Phosphoryl, ...). Layer-2 Subtype is still in the
        # records and shown in the per-ligand detail view; here we filter at the broad level
        # to keep the dropdown manageable. See vocab/chemical_families.json for the taxonomy.
        family_options = sorted({
            getattr(l, 'chemical_family', None) for l in db.get_all_ligands()
            if getattr(l, 'chemical_family', None)
        })
        extractant_types_filter = st.multiselect(
            "Ligand Subtype",
            options=family_options,
            default=[],
            help="Filter by Subtype (21 categories: DGA, Malonamide, DAPhen, Phosphine_oxide, ...). See vocab/chemical_families.json for the taxonomy."
        )

        # Mechanism filter - HIDDEN (data quality concerns)
        # mechanisms = st.multiselect(
        #     "Extraction Mechanism",
        #     options=["Solvation", "Cation_Exchange"],
        #     default=[],
        #     help="Filter by extraction mechanism"
        # )
        mechanisms = []  # Disabled

        # Source filters exist only where records carry source labels (not in the public release)
        has_sources = db.has_source_labels()

        # Liu data only toggle
        show_liu_only = has_sources and st.checkbox(
            "📊 Show only JACSAu-2022 data",
            value=False,
            help="Filter to show only ligands from Liu et al. 2022 JACS Au dataset"
        )

        # Trend filter
        trends = st.multiselect(
            "Selectivity Trend",
            options=list(TREND_DISPLAY_ORDER[:-1]),  # every value except "Unknown"
            default=[],
            help="Filter by selectivity trend classification"
        )

        # Data source filter
        sources = st.multiselect(
            "Data Source",
            options=["Liu_2022_JACSAu", "SAFE_01262026", "OECD_NEA_IDEAL", "Pramanik_2024_JACS",
                     "Literature_curated"],
            default=[],
            help="Filter by data source (can select multiple)"
        ) if has_sources else []

        # Element filter
        with st.expander("Required Elements"):
            required_elements = st.multiselect(
                "Must have data for:",
                LANTHANIDES,
                default=[],
                help="Filter ligands that have data for specific elements"
            )

        # LogD Range Filter
        with st.expander("📊 logD Range Filter"):
            logd_range = st.slider(
                "Select log(D) range",
                min_value=-5.0,
                max_value=7.0,
                value=(-4.0, 3.0),
                step=0.1,
                help=(
                    "Filter ligands by log(Distribution Ratio) range.\n"
                    "logD = log(concentration_organic / concentration_aqueous)\n\n"
                    "• logD < 0: Poor extraction (prefers aqueous phase)\n"
                    "• logD = 0: Equal distribution\n"
                    "• logD > 0: Strong extraction (prefers organic phase)\n"
                    "• logD > 2: Excellent extraction (>99% to organic)"
                )
            )
            logd_min, logd_max = logd_range

            # Show filter status
            if logd_min != -4.0 or logd_max != 3.0:
                st.caption(f"✓ Filtering: {logd_min:.1f} ≤ logD ≤ {logd_max:.1f}")
                if logd_min == -5.0 and logd_max == 7.0:
                    st.caption("ℹ️ Full range selected (all data)")
            else:
                st.caption("ℹ️ Default range (98% of data)")

        st.markdown("---")

        #=====================================================================
        # RESULTS COUNT BADGE & JUMP BUTTON
        #=====================================================================

        # Calculate filtered ligands count (duplicate filter logic for real-time feedback)
        ligands_for_count = db.get_all_ligands()

        if search_text:
            ligands_for_count = [l for l in ligands_for_count if search_text.lower() in l.name.lower()]

        if extractant_types_filter:
            ligands_for_count = [l for l in ligands_for_count if getattr(l, 'chemical_family', None) in extractant_types_filter]

        if mechanisms:
            ligands_for_count = [l for l in ligands_for_count if l.mechanism in mechanisms]

        if trends:
            ligands_for_count = [l for l in ligands_for_count if l.selectivity_trend in trends]

        if show_liu_only:
            ligands_for_count = [l for l in ligands_for_count if 'Liu_2022_JACSAu' in l.data_sources]

        if sources and not show_liu_only:
            ligands_for_count = [l for l in ligands_for_count if any(src in l.data_sources for src in sources)]

        if required_elements:
            ligands_for_count = [l for l in ligands_for_count if all(
                any(elem in c.get('distribution_ratios', {}) for c in l.experimental_conditions)
                for elem in required_elements
            )]

        # Apply logD range filter
        if logd_min != -4.0 or logd_max != 3.0:
            ligands_for_count = [lig for lig in ligands_for_count
                                 if has_value_in_range(lig, logd_min, logd_max)]

        ligand_count = len(ligands_for_count)

        # Color-coded results count
        if ligand_count >= 50:
            count_color = "🟢"
        elif ligand_count >= 10:
            count_color = "🟡"
        else:
            count_color = "🔴"

        st.markdown(f"### {count_color} **{ligand_count} ligands found**")

        # Active filter indicator
        active_filters = []
        if search_text:
            active_filters.append("Search")
        if extractant_types_filter:
            active_filters.append("Type")
        if mechanisms:
            active_filters.append("Mechanism")
        if trends:
            active_filters.append("Selectivity")
        if show_liu_only:
            active_filters.append("JACSAu-2022 Only")
        if sources:
            active_filters.append("Data Source")
        if logd_min != -4.0 or logd_max != 3.0:
            active_filters.append("logD Range")
        if required_elements:
            active_filters.append(f"Required Elements ({len(required_elements)})")

        if active_filters:
            st.caption(f"✓ **{len(active_filters)} filters active:** {', '.join(active_filters)}")
        else:
            st.caption("ℹ️ No filters applied (showing all ligands)")

        # Jump to Results button with anchor link
        st.markdown(
            f"""
            <style>
            .jump-button {{
                display: block;
                text-align: center;
                padding: 14px 20px;
                margin: 10px 0;
                background: linear-gradient(135deg, #0068c9 0%, #0056b3 100%);
                color: white;
                border: 2px solid rgba(255, 255, 255, 0.3);
                border-radius: 8px;
                text-decoration: none;
                font-weight: 700;
                font-size: 16px;
                letter-spacing: 0.5px;
                text-shadow: 0 1px 2px rgba(0, 0, 0, 0.2);
                box-shadow: 0 2px 4px rgba(0, 0, 0, 0.15);
                transition: all 0.3s ease;
            }}
            .jump-button:hover {{
                transform: translateY(-2px);
                background: linear-gradient(135deg, #1a7dd9 0%, #1566c3 100%);
                border-color: rgba(255, 255, 255, 0.5);
                box-shadow: 0 4px 12px rgba(0, 0, 0, 0.25);
                text-decoration: none;
                color: white;
            }}
            </style>
            <a href="#ligand-browser" class="jump-button">
                ↓ View Filtered Ligands ({ligand_count})
            </a>
            """,
            unsafe_allow_html=True
        )

        st.markdown("---")

        # Database statistics
        with st.expander("📊 Database Statistics"):
            stats = db.get_statistics()
            all_ligands = db.get_all_ligands()

            # Calculate total conditions and merged ligands
            total_conditions = sum(len(lig.experimental_conditions) for lig in all_ligands)
            merged_ligands = sum(
                1 for lig in all_ligands
                if hasattr(lig, 'metadata') and isinstance(lig.metadata, dict)
                and isinstance(lig.metadata.get('data_sources', []), list)
                and len(lig.metadata.get('data_sources', [])) > 1
            )

            st.metric("Total Ligands", stats['total_ligands'])
            st.metric("Total Conditions", total_conditions)
            if db.has_source_labels():
                st.metric("Merged Ligands", merged_ligands)

            # Ligand family breakdown (Layer-1, 13 broad categories)
            st.write("**By Ligand Family:**")
            type_counts = {}
            for lig in all_ligands:
                type_name = getattr(lig, 'chemical_family', None) or lig.extractant_type or 'Unknown'
                type_counts[type_name] = type_counts.get(type_name, 0) + 1

            for etype in sorted(type_counts.keys(), key=lambda x: type_counts[x], reverse=True):
                count = type_counts[etype]
                pct = 100 * count / len(all_ligands)
                st.write(f"- {etype}: {count} ({pct:.1f}%)")

            st.write("**By Trend:**")
            for trend, count in stats.get('trends', {}).items():
                st.write(f"- {trend}: {count}")

        # Display Settings
        st.markdown("---")
        st.markdown("### 🖥️ Display Settings")

        # Initialize session state for fullscreen fonts
        if 'fullscreen_fonts' not in st.session_state:
            st.session_state.fullscreen_fonts = False

        # Toggle button
        st.session_state.fullscreen_fonts = st.checkbox(
            "Full Screen Fonts",
            value=st.session_state.fullscreen_fonts,
            help="Enable larger font sizes optimized for full-screen viewing"
        )

    return (page, show_liu_only, search_text, extractant_types_filter,
            mechanisms, trends, sources, required_elements, logd_min, logd_max)


def _render_overview_tab(db, show_liu_only):
    """
    Render the Database Overview, Lanthanide Coverage and Eu/La sections.

    Owns the ligand list because the Liu-only filter emits a status message
    here, so the filtering has to happen in render order.

    Returns:
        list: the (possibly Liu-filtered) ligands, reused by the other tabs
    """
    st.markdown("---")
    st.subheader("📊 Database Overview")

    # Get all ligands for overview statistics
    all_ligands = db.get_all_ligands()

    # Apply Liu-only filter to all visualizations if checked
    if show_liu_only:
        all_ligands = [l for l in all_ligands if 'Liu_2022_JACSAu' in l.data_sources]
        st.info(f"📊 Showing JACSAu-2022 data only: {len(all_ligands)} ligands")

    # The subtype bar takes full width (21 categories need the horizontal room).
    # Selectivity-trend pie sits below it on its own row, slightly inset.
    fig_type = plot_ligand_type_pie(all_ligands)
    st.plotly_chart(fig_type, use_container_width=True)

    trend_left, trend_mid, trend_right = st.columns([1, 2, 1])
    with trend_mid:
        fig_trend = plot_selectivity_trend_pie(all_ligands)
        st.plotly_chart(fig_trend, use_container_width=True)

    # Lanthanide Distribution Donut Chart
    st.markdown("---")
    st.subheader("🔬 Lanthanide Coverage")

    col3, col4 = st.columns([1, 1])

    with col3:
        fig_lanthanides = plot_lanthanide_distribution_donut(all_ligands)
        st.plotly_chart(fig_lanthanides, use_container_width=True)
        st.caption("Shows how many ligands have measurements for each lanthanide element.")

    with col4:
        # Show top 5 most measured and least measured
        element_counts = {}
        for lig in all_ligands:
            elements_in_ligand = set()
            for cond in lig.experimental_conditions:
                dist_ratios = cond.get('distribution_ratios', {})
                for elem in LANTHANIDES:
                    if elem in dist_ratios:
                        elements_in_ligand.add(elem)
            for elem in elements_in_ligand:
                element_counts[elem] = element_counts.get(elem, 0) + 1

        if element_counts:
            sorted_elements = sorted(element_counts.items(), key=lambda x: x[1], reverse=True)

            st.write("**Most Measured Lanthanides:**")
            for elem, count in sorted_elements[:5]:
                percentage = (count / len(all_ligands)) * 100
                st.write(f"- **{elem}**: {count} ligands ({percentage:.1f}%)")

            st.write("")
            st.write("**Least Measured Lanthanides:**")
            for elem, count in sorted_elements[-5:]:
                percentage = (count / len(all_ligands)) * 100
                st.write(f"- **{elem}**: {count} ligands ({percentage:.1f}%)")

    # Eu/La Separation Factor Scatter Plots
    st.markdown("---")
    st.subheader("⚛️ Eu/La Separation Analysis")

    fig_eu_la = plot_eu_la_scatter(all_ligands)
    st.plotly_chart(fig_eu_la, use_container_width=True)

    # Add explanation
    with st.expander("ℹ️ About Eu/La Separation Factors"):
        st.write("""
        **Separation Factor Definition:**
        - log(SF[La/Eu]) = log(D[La]) - log(D[Eu])
        - **Positive SF**: La-selective (prefers Lanthanum)
        - **Negative SF**: Eu-selective (prefers Europium)
        - **SF ≈ 0**: No selectivity between Eu and La

        **Why Eu and La?**
        - Europium (Eu) and Lanthanum (La) are the most frequently measured lanthanides in this database
        - They represent light (La) and middle (Eu) rare earth elements
        - Separation of these elements is important for rare earth processing

        **Interpreting the Plots:**
        - Left plot: Shows how Eu extraction correlates with selectivity
        - Right plot: Shows how La extraction correlates with selectivity
        - Points above the dashed line favor La extraction
        - Points below the dashed line favor Eu extraction
        """)

    st.markdown("---")

    return all_ligands


def _render_single_element_tab(db, all_ligands):
    """Render the Single Element Distribution Viewer."""
    #=========================================================================
    # Single Element logD Viewer
    #=========================================================================

    # Initialize session state for selected element
    if 'selected_element' not in st.session_state:
        st.session_state.selected_element = 'Eu'  # Default to Europium

    # Custom CSS for periodic table element buttons
    st.markdown("""
    <style>
        /* Periodic table element button styling */
        div[data-testid="column"] button {
            width: 100% !important;
            padding: 8px 4px !important;
            font-size: 16px !important;
            font-weight: bold !important;
            border: 2px solid #333 !important;
            border-radius: 4px !important;
            background: linear-gradient(135deg, #B8E6F0 0%, #A0D8E8 100%) !important;
            color: #333 !important;
            transition: all 0.2s ease !important;
            height: 60px !important;
            display: flex !important;
            flex-direction: column !important;
            align-items: center !important;
            justify-content: center !important;
        }

        div[data-testid="column"] button:hover {
            transform: scale(1.05) !important;
            border-color: #0066CC !important;
            box-shadow: 0 4px 8px rgba(0,0,0,0.2) !important;
        }

        /* Selected element styling */
        div[data-testid="column"] button[kind="primary"] {
            background: linear-gradient(135deg, #4A90E2 0%, #357ABD 100%) !important;
            color: white !important;
            border-color: #0066CC !important;
            box-shadow: 0 0 12px rgba(74, 144, 226, 0.6) !important;
        }

        /* Atomic number styling */
        .element-atomic-number {
            font-size: 10px;
            opacity: 0.7;
        }
    </style>
    """, unsafe_allow_html=True)

    st.subheader("🔬 Single Element Distribution Viewer")

    # Periodic Table Element Selector
    selected_element = render_element_selector(
        state_key='selected_element',
        button_prefix='elem',
        label="**Select Lanthanide Element:**",
        caption_prefix="Selected: ",
    )

    # X-Axis selector remains below the periodic table
    st.markdown("**Visualization Options:**")
    x_axis_option = st.radio(
        "X-Axis",
        options=["Ligand Types", "Solvents", "Solvent Categories", "Ligand Acidity", "Acid Concentration", "Extractant Concentration", "Temperature"],
        index=0,
        horizontal=True,
        help="Choose what to plot on the x-axis"
    )

    # Show category info when Solvent Categories is selected
    if x_axis_option == "Solvent Categories":
        with st.expander("ℹ️ Solvent Category Definitions"):
            st.markdown("""
            **Solvent categories group similar organic diluents:**
            - **aliphatic**: Non-polar saturated hydrocarbons (n-dodecane, Isopar L)
            - **aromatic**: Moderate polarity benzene-based (toluene, DIPB, tBuB)
            - **chlorinated_hydrocarbon**: High polarity halogenated (chloroform, DCM, CCl4)
            - **alcohol**: Phase modifiers with hydroxyl groups (1-octanol, Isodecanol)
            - **petroleum_mixture**: Industrial blends (kerosene, sulfonated kerosene)
            - **organophosphorus**: Dual extractant/diluent (TBP)
            - **ketone**: Polar aprotic carbonyl-containing (cyclohexanone)
            - **aromatic_nitro**: High polarity nitro-substituted (nitrobenzene)
            - **aromatic_fluoro**: Fluorinated aromatics (meta-nitrobenzotrifluoride)
            - **aromatic_sulfone**: Sulfone-containing (phenyl trifluoromethyl sulfone)
            - **Other/Unknown**: Unmapped or unknown solvents

            Categories aggregate measurements across multiple individual solvents for easier pattern recognition.
            """)

    # Show acidity info when Ligand Acidity is selected
    if x_axis_option == "Ligand Acidity":
        with st.expander("ℹ️ Ligand Acidity Classification"):
            st.markdown("""
            **Ligand acidity refers to the extraction mechanism:**

            **Acidic Ligands** (Red) - 4 ligands (2.3%):
            - Release H⁺ during metal extraction (cation exchange mechanism)
            - M³⁺(aq) + 3HL(org) ⇌ ML₃(org) + 3H⁺(aq)
            - Examples:
              - **Carboxylic acids**: HEDTA (with -COOH groups)
              - **Phosphonic/phosphoric acids**: HDHP (with P(=O)-OH groups)
              - **Sulfonic acids**: rare in this database

            **Neutral Ligands** (Blue) - 168 ligands (94.9%):
            - Coordinate metal without H⁺ release (solvation mechanism)
            - M³⁺(aq) + 3A⁻(aq) + nL(org) ⇌ MA₃·nL(org)
            - Examples:
              - **Diglycolamides**: TODGA, TEHDGA, THDGA (DGA family)
              - **Malonamides**: DMDOHEMA, DMDBTDMA
              - **Phosphine oxides**: TOPO, CMPO, TBPO (R-P=O)
              - **Phosphate esters**: POPPhen derivatives (P(=O)-O-R)
              - **Aromatic heterocycles**: BTP, phenanthroline
              - **Crown ethers**: 18-crown-6, benzo-15-crown-5
              - **Amides**: DOHyA (with -CONH- groups)

            **Unknown** (Gray) - 5 ligands (2.8%):
            - Missing SMILES structure data (cannot classify without molecular structure)

            **Classification Method:**
            1. **SMARTS pattern matching** (primary, 97.2% success rate):
               - Detects functional groups from molecular structure (SMILES)
               - Acidic: Requires -OH with H (e.g., -COOH, P(=O)-OH)
               - Neutral: Amides, phosphine oxides, phosphate esters, heterocycles
            2. **Extractant type** (secondary):
               - Based on structure-validated classifications (97.3% accuracy)
            3. **Keyword matching** (fallback):
               - Checks ligand names for known acidic/neutral keywords

            **Recent Improvement**: SMARTS-based classification increased accuracy from 93.2% → 97.2%
            """)

    # Extract data for selected element
    element_df = extract_element_data(db, selected_element, all_ligands)

    # Solvent filter (only show if we have data)
    if not element_df.empty:
        # Get unique solvents from the data
        unique_solvents = sorted(element_df['solvent'].unique())

        # Solvent multiselect
        selected_solvents = st.multiselect(
            "Filter by Solvent",
            options=unique_solvents,
            default=[],
            help="Select one or more solvents to filter the data. Leave empty to show all solvents.",
            placeholder="All solvents"
        )

        # Apply solvent filter
        if selected_solvents:
            element_df = element_df[element_df['solvent'].isin(selected_solvents)]

            if element_df.empty:
                st.warning(f"No data for {selected_element} with the selected solvent(s). Showing all data instead.")
                # Reset to show all
                element_df = extract_element_data(db, selected_element, all_ligands)

    # Create and display plot
    if not element_df.empty:
        fig_element = plot_element_distribution(element_df, selected_element, x_axis_option)
        st.plotly_chart(fig_element, use_container_width=True)
        n_bounds = int(element_df['censor'].notna().sum())
        if n_bounds:
            st.caption(f"▲ D > x and ▼ D < x: {n_bounds} bound(s), where the paper gives only a limit. "
                       "They are drawn at the limit but left out of the box, violin and statistics.")

        # Explanation
        with st.expander("ℹ️ How to Use This Viewer"):
            st.write(f"""
            **About this visualization:**
            - Shows all measured distribution coefficients (D-values) for **{selected_element}** in the database
            - Distribution coefficient: D = [M]org / [M]aq (concentration in organic phase / aqueous phase)
            - Higher log(D) values indicate better extraction into the organic phase

            **X-Axis Options:**
            - **Ligand Types**: Box plot showing distribution of log(D) across ligand families (DGA, Phosphine oxide, etc.)
            - **Solvents**: Box plot showing distribution of log(D) across organic solvents (~80 unique)
            - **Solvent Categories**: Box plot showing distribution across solvent chemical types (aliphatic, aromatic, chlorinated, etc.)
            - **Ligand Acidity**: Box plot comparing acidic vs neutral ligands (cation exchange vs solvation mechanisms)
            - **Acid Concentration**: Scatter plot showing relationship with acidity
            - **Extractant Concentration**: Scatter plot showing loading effects
            - **Temperature**: Scatter plot showing thermodynamic trends

            **Box + Violin Plot Features (for Ligand Types, Solvents, Solvent Categories, Ligand Acidity):**
            - **Violin plot** (background): Semi-transparent colored shape showing probability density distribution
              - Wider regions = more data points at that log(D) value
              - Narrow regions = fewer data points
              - Provides full picture of data distribution
            - **Box plot** (foreground): Overlaid on violin for statistical summary
              - **Box**: Shows Q1 (25th percentile) to Q3 (75th percentile)
              - **Median line**: Thick line inside box (50th percentile)
              - **Mean line**: Thin line inside box (average value)
              - **Whiskers**: Extend to min/max within 1.5×IQR
              - **Outliers**: Individual colored points beyond whiskers
            - **Color coding**: Each category has a distinct color for easy visual differentiation
              - Ligand Types: Blue (Diglycolamide), Purple (Aromatic_Heterocycle), Orange (Organophosphorus), etc.
              - Solvents: Blue (n-dodecane), Orange (toluene), Green (chloroform), etc.
              - Solvent Categories: Blue (aliphatic), Orange (aromatic), Green (chlorinated), etc.
              - Ligand Acidity: Red (Acidic - cation exchange), Blue (Neutral - solvation)
            - **Hover**: Shows all statistics (median, mean, quartiles, count, top ligands)
            - Hover over outlier points to see individual measurement details (ligand name, solvent, conditions)

            **Solvent Categories:**
            Solvent categories group similar organic diluents by chemical type:
            - **aliphatic**: Non-polar saturated hydrocarbons (n-dodecane, Isopar L)
            - **aromatic**: Moderate polarity benzene-based (toluene, benzene, DIPB, tBuB)
            - **chlorinated_hydrocarbon**: High polarity halogenated (chloroform, DCM, CCl4)
            - **alcohol**: Phase modifiers with hydroxyl groups (1-octanol, Isodecanol)
            - **petroleum_mixture**: Industrial blends (kerosene, sulfonated kerosene)
            - **organophosphorus**: Dual extractant/diluent (TBP)
            - **ketone**: Polar aprotic carbonyl-containing (cyclohexanone)
            - **aromatic_nitro**: High polarity nitro-substituted (nitrobenzene)
            - **aromatic_fluoro**: Fluorinated aromatics (meta-nitrobenzotrifluoride)
            - **aromatic_sulfone**: Sulfone-containing (phenyl trifluoromethyl sulfone)
            - **Other/Unknown**: Unmapped or unknown solvents

            **Solvent Filter:**
            - Filter data by specific organic solvent(s)
            - JACSAu-2022 data now has 21+ specific solvents (n-dodecane, toluene, chloroform, etc.)
            - SAFE and IDEAL data also have solvent information
            - Use this to compare extractant performance across different diluents

            **Tips:**
            - Hover over box plots to see comprehensive statistics (median, mean, quartiles, range, top ligands)
            - Hover over outlier points to see individual measurement details
            - Hover over scatter points to see ligand name, solvent, and experimental conditions
            - Box plots are sorted left-to-right by median log(D) (best performing categories on the left)
            - Use this to identify the best extractant types for {selected_element}
            - Compare with other elements by changing the selection above
            - Try filtering by common solvents like "n-dodecane" or "toluene" to see solvent effects
            - Use "Solvent Categories" for a cleaner view when many solvents are present
            """)
    else:
        st.warning(f"No measurements found for {selected_element} in the database.")
        st.info("Try selecting a different element from the dropdown above.")

    st.markdown("---")


def _render_logsf_tab(db, all_ligands):
    """Render the logSF (Separation Factor) Viewer."""
    #=========================================================================
    # logSF (Separation Factor) Viewer
    #=========================================================================

    # Initialize session state for element selections
    if 'selected_element_ln1' not in st.session_state:
        st.session_state.selected_element_ln1 = 'Eu'
    if 'selected_element_ln2' not in st.session_state:
        st.session_state.selected_element_ln2 = 'La'

    st.subheader("⚖️ logSF (Separation Factor) Viewer")

    # Explanation
    st.markdown("""
    Visualize separation factors between two lanthanides.
    **logSF[Ln1/Ln2] = logD[Ln1] - logD[Ln2]**
    - **Positive logSF**: Ln1 is preferentially extracted
    - **Negative logSF**: Ln2 is preferentially extracted
    - **logSF ≈ 0**: No selectivity
    """)

    # First Element Selector (Ln1)
    render_element_selector(
        state_key='selected_element_ln1',
        button_prefix='ln1',
        label="**Select First Element (Ln1):**",
        caption_prefix="Ln1: ",
    )

    # Second Element Selector (Ln2)
    render_element_selector(
        state_key='selected_element_ln2',
        button_prefix='ln2',
        label="**Select Second Element (Ln2):**",
        caption_prefix="Ln2: ",
    )

    # Store in variables
    ln1 = st.session_state.selected_element_ln1
    ln2 = st.session_state.selected_element_ln2

    # Validate selection
    if ln1 == ln2:
        st.warning(f"⚠️ Please select two different elements. Currently both are set to {ln1}.")
    else:
        st.info(f"**Current Separation:** logSF[{ln1}/{ln2}] = logD[{ln1}] - logD[{ln2}]")

        # X-Axis selector
        st.markdown("**Visualization Options:**")
        x_axis_option_sf = st.radio(
            "X-Axis",
            options=["Ligand Types", "Solvents", "Solvent Categories",
                     "Ligand Acidity", "Acid Concentration",
                     "Extractant Concentration", "Temperature"],
            index=0,
            horizontal=True,
            key="x_axis_sf"
        )

        # Category info expanders (when relevant)
        if x_axis_option_sf == "Solvent Categories":
            with st.expander("ℹ️ Solvent Category Definitions"):
                st.markdown("""
                **Solvent categories group similar organic diluents:**
                - **aliphatic**: Non-polar saturated hydrocarbons (n-dodecane, Isopar L)
                - **aromatic**: Moderate polarity benzene-based (toluene, DIPB, tBuB)
                - **chlorinated_hydrocarbon**: High polarity halogenated (chloroform, DCM, CCl4)
                - **alcohol**: Phase modifiers with hydroxyl groups (1-octanol, Isodecanol)
                - **petroleum_mixture**: Industrial blends (kerosene, sulfonated kerosene)
                - **organophosphorus**: Dual extractant/diluent (TBP)
                - **ketone**: Polar aprotic carbonyl-containing (cyclohexanone)
                - **aromatic_nitro**: High polarity nitro-substituted (nitrobenzene)
                - **aromatic_fluoro**: Fluorinated aromatics (meta-nitrobenzotrifluoride)
                - **aromatic_sulfone**: Sulfone-containing (phenyl trifluoromethyl sulfone)
                - **Other/Unknown**: Unmapped or unknown solvents

                Categories aggregate measurements across multiple individual solvents for easier pattern recognition.
                """)

        if x_axis_option_sf == "Ligand Acidity":
            with st.expander("ℹ️ Ligand Acidity Classification"):
                st.markdown("""
                **Ligand acidity refers to the extraction mechanism:**

                **Acidic Ligands** (Red):
                - Release H⁺ during metal extraction (cation exchange mechanism)
                - M³⁺(aq) + 3HL(org) ⇌ ML₃(org) + 3H⁺(aq)
                - Examples: DEHPA, D2EHPA, HDEHP, HEDTA, versatic acid
                - Carboxylic acids, phosphonic/phosphoric acids, sulfonic acids

                **Neutral Ligands** (Blue):
                - Coordinate without H⁺ release (solvation mechanism)
                - M³⁺(aq) + 3L(org) ⇌ ML₃(org)
                - Examples: DGA family (TODGA, TEHDGA), malonamides, phosphine oxides (TOPO, CMPO)
                - Most ligands in database (diglycolamides, aromatic heterocycles)

                **Unknown** (Gray):
                - Insufficient structure/mechanism data for classification
                """)

        # Extract data
        sf_df = extract_separation_factor_data(db, ln1, ln2, all_ligands)
        sf_no_information = sf_df.attrs.get('no_information', 0)

        # Solvent Filter
        if not sf_df.empty:
            st.markdown("**Filters:**")
            available_solvents_sf = sorted(sf_df['solvent'].unique())

            with st.expander("🧪 Filter by Solvent"):
                selected_solvents_sf = st.multiselect(
                    "Select solvents to include:",
                    available_solvents_sf,
                    default=available_solvents_sf,
                    key="solvent_filter_sf"
                )

                # Apply filter
                if selected_solvents_sf and selected_solvents_sf != available_solvents_sf:
                    sf_df = sf_df[sf_df['solvent'].isin(selected_solvents_sf)]

        # Create and display plot
        if not sf_df.empty:
            fig_sf = plot_separation_factor_distribution(sf_df, ln1, ln2, x_axis_option_sf)
            st.plotly_chart(fig_sf, use_container_width=True)
            n_sf_bounds = int(sf_df['censor'].notna().sum())
            if n_sf_bounds or sf_no_information:
                st.caption(
                    (f"▲ logSF > x and ▼ logSF < x: {n_sf_bounds} pair(s) where one D is only a limit, "
                     "drawn at the limit but left out of the box, violin and statistics. " if n_sf_bounds else "")
                    + (f"{sf_no_information} pair(s) with both D's limited the same way say nothing about "
                       "the ratio and are not shown." if sf_no_information else ""))

            # Documentation
            with st.expander("ℹ️ How to Use This Viewer"):
                st.write(f"""
                **About this visualization:**
                - Shows separation factors between **{ln1}** and **{ln2}**
                - Only includes conditions where BOTH elements were measured
                - Positive values indicate {ln1} preference (logD[{ln1}] > logD[{ln2}])
                - Negative values indicate {ln2} preference (logD[{ln2}] > logD[{ln1}])
                - Gray dashed line at y=0 shows no selectivity

                **X-Axis Options:**
                - **Ligand Types**: Compare extractant families (DGA, phosphine oxides, malonamides, etc.)
                - **Solvents**: View performance across different organic diluents
                - **Solvent Categories**: Grouped by solvent class (aliphatic, aromatic, chlorinated, etc.)
                - **Ligand Acidity**: Compare acidic (cation exchange) vs neutral (solvation) ligands
                - **Acid Concentration**: How separation changes with aqueous acid molarity
                - **Extractant Concentration**: Effect of organic phase extractant concentration
                - **Temperature**: Temperature dependence of separation

                **Box + Violin Plot Features:**
                - **Violin (colored background)**: Probability density showing data distribution shape
                - **Box (white overlay)**: Statistical summary (median, quartiles, outliers)
                - **Hover over box**: See statistics and top-performing ligands for that category
                - **Hover over points**: See individual measurement details (ligand name, conditions, logD values)

                **Tips:**
                - Box plots are sorted left-to-right by median logSF (best performing categories on the left)
                - Use this to identify the best extractants for {ln1}/{ln2} separation
                - Try different element pairs to explore selectivity patterns
                - Use solvent filter to examine solvent effects
                - Look for categories with high positive logSF for {ln1}-selective ligands
                - Look for categories with high negative logSF for {ln2}-selective ligands
                """)
        else:
            st.warning(f"⚠️ No paired measurements found for {ln1} and {ln2} under the current filter settings.")
            st.info("💡 Try selecting different elements or adjusting the solvent filter to see more data.")

    st.markdown("---")


def _filter_ligands(db, search_text, extractant_types_filter, mechanisms, trends,
                    show_liu_only, sources, required_elements, logd_min, logd_max):
    """
    Apply every sidebar filter to the ligand list.

    Pure logic with no Streamlit calls, so it can be unit-tested directly.

    Returns:
        list: the ligands passing all active filters
    """
    # Apply filters
    ligands = db.get_all_ligands()

    if search_text:
        ligands = [l for l in ligands if search_text.lower() in l.name.lower()]

    if extractant_types_filter:
        ligands = [l for l in ligands if getattr(l, 'chemical_family', None) in extractant_types_filter]

    if mechanisms:
        ligands = [l for l in ligands if l.mechanism in mechanisms]

    if trends:
        ligands = [l for l in ligands if l.selectivity_trend in trends]

    # Apply Liu-only filter
    if show_liu_only:
        ligands = [l for l in ligands if 'Liu_2022_JACSAu' in l.data_sources]

    # Apply data source filter (if not using Liu-only)
    if sources and not show_liu_only:
        # Filter ligands that have at least one of the selected sources
        ligands = [l for l in ligands if any(src in l.data_sources for src in sources)]

    if required_elements:
        # Filter ligands that have data for all required elements
        ligands = [l for l in ligands if all(
            any(elem in c.get('distribution_ratios', {}) for c in l.experimental_conditions)
            for elem in required_elements
        )]

    # Apply logD range filter
    if logd_min is not None and logd_max is not None:
        # Only filter if range is not the full default (-4.0 to 3.0)
        # This avoids unnecessary computation when user hasn't adjusted slider
        if logd_min != -4.0 or logd_max != 3.0:
            # at least one measured value (not a bound) in the logD range
            ligands = [lig for lig in ligands if has_value_in_range(lig, logd_min, logd_max)]

    return ligands


def _render_ligand_table(ligands):
    """
    Render the browser table and its CSV download button.

    Returns:
        The st.dataframe selection event; the selected row drives the detail view.
    """
    df_data = []
    for lig in ligands:
        dq = lig.data_quality or {}
        quality_str = f"{dq.get('overall_quality_score', '?')} ({dq.get('quality_grade', '?')})" if dq else 'N/A'
        df_data.append({
            'Name': lig.name,
            'SMILES': lig.smiles[:30] + '...' if len(lig.smiles or '') > 30 else (lig.smiles or 'N/A'),
            # 'Mechanism': lig.mechanism,  # Hidden due to data quality concerns
            'Ligand Family': getattr(lig, 'chemical_family', None) or lig.extractant_type or 'N/A',
            'Subtype': getattr(lig, 'subtype', None) or lig.community_name or 'N/A',
            'Trend': lig.selectivity_trend,
            'Confidence': lig.trend_confidence or 'N/A',
            'Quality': quality_str,
            'Conditions': len(lig.experimental_conditions),
            'DOI': _get_primary_doi(lig) or '—',
        })

    df = pd.DataFrame(df_data)

    # Anchor point for jump navigation
    st.markdown('<div id="ligand-browser"></div>', unsafe_allow_html=True)

    # Header with back-to-filters link
    col1, col2 = st.columns([4, 1])
    with col1:
        st.subheader(f"Ligand Browser ({len(ligands)} ligands found)")
    with col2:
        st.markdown(
            '<a href="#filters" style="'
            'display:block; text-align:right; padding:8px 0; '
            'color:#0068c9; text-decoration:none; font-size:14px; '
            'font-weight:500;" '
            'title="Scroll back to filter controls">'
            '↑ Refine Filters'
            '</a>',
            unsafe_allow_html=True
        )

    # Display table with selection
    event = st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        on_select="rerun",
        selection_mode="single-row",
        key="ligand_table"
    )

    return event


def _render_ligand_identity(selected_ligand):
    """Left detail column: structure, chemical identifiers and references."""
    st.subheader("Structure")

    if selected_ligand.smiles and selected_ligand.smiles != '-':
        img = render_molecule(selected_ligand.smiles)
        if img:
            st.image(img, use_container_width=True)
        elif not RDKIT_AVAILABLE:
            st.info(
                "Structure drawing needs RDKit, which is not installed in this "
                "deployment. The SMILES below describes the same molecule."
            )
        else:
            st.info("This SMILES could not be rendered.")

        st.caption(f"**SMILES**: `{selected_ligand.smiles}`")
    else:
        st.warning("No valid SMILES available for this ligand")

    # Additional info
    # st.caption(f"**Mechanism**: {selected_ligand.mechanism}")  # Hidden due to data quality concerns
    _family = getattr(selected_ligand, 'chemical_family', None) or selected_ligand.extractant_type or 'N/A'
    _subtype = getattr(selected_ligand, 'subtype', None) or getattr(selected_ligand, 'community_name', None) or 'N/A'
    st.caption(f"**Ligand family**: {_family}  •  **Subtype**: {_subtype}")

    if hasattr(selected_ligand, 'alternative_names') and selected_ligand.alternative_names:
        st.caption(f"**Also known as**: {', '.join(selected_ligand.alternative_names)}")

    # Chemical identifiers (added by IUPAC-naming pipeline, v8_named)
    if getattr(selected_ligand, 'iupac_name', None):
        st.divider()
        st.markdown("**Chemical Identifiers**")
        # The public release drops the naming provenance (Name_Source, Name_Confidence):
        # show the badge only where the record carries it, never "unknown, unknown".
        if selected_ligand.name_source or selected_ligand.name_confidence:
            conf_emoji = {"high": "🟢", "medium": "🟡", "low": "🟠"}.get(
                selected_ligand.name_confidence or "", "⚪"
            )
            st.markdown(
                f"**IUPAC name** {conf_emoji} "
                f"<span style='font-size:0.85em'>({selected_ligand.name_source or 'unknown'}, "
                f"{selected_ligand.name_confidence or 'unknown'} confidence)</span>",
                unsafe_allow_html=True,
            )
        else:
            st.markdown("**IUPAC name**")
        st.code(selected_ligand.iupac_name, language=None)
        if selected_ligand.pubchem_cid:
            cid_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{selected_ligand.pubchem_cid}"
            st.caption(f"**PubChem CID**: [{selected_ligand.pubchem_cid}]({cid_url})")
        if selected_ligand.cas_number:
            st.caption(f"**CAS**: `{selected_ligand.cas_number}`")
        if selected_ligand.inchi_key:
            st.caption(f"**InChI Key**: `{selected_ligand.inchi_key}`")

    # Data Sources and References
    st.divider()
    st.subheader("Data Sources & References")

    # Get references from metadata (new structure)
    if hasattr(selected_ligand, 'metadata') and 'references' in selected_ligand.metadata:
        metadata_refs = selected_ligand.metadata['references']

        for ref in metadata_refs:
            # Display source database
            source_db = ref.get('source_database', 'Unknown')
            st.write(f"**{source_db}**")

            # Display specific reference info for Liu 2022
            if 'reference_number' in ref:
                ref_num = ref['reference_number']
                st.markdown(f"📑 Reference {ref_num}")

            # Display DOI if available
            if ref.get('doi'):
                st.markdown(f"📄 DOI: [{ref['doi']}]({doi_url(ref['doi'])})")

            # Display citation
            if ref.get('citation'):
                st.caption(f"*{ref['citation']}*")

            # Display URL for online databases
            if ref.get('url') and 'doi' not in ref.get('url', ''):
                st.markdown(f"🔗 [Database Link]({ref['url']})")

            # Display notes
            if ref.get('notes'):
                st.caption(f"ℹ️ {ref['notes']}")

            st.write("")  # Spacing

    # Per-condition provenance (schema v0.3): one entry per measuring paper.
    else:
        papers = {}
        for cond in selected_ligand.experimental_conditions:
            ref = cond.get('reference') or {}
            doi = primary_doi(ref)
            key = doi or ref.get('URL')
            if not key:
                continue
            entry = papers.setdefault(key, {
                'doi': doi, 'url': ref.get('URL'), 'n': 0,
                'locators': [], 'aggregators': [], 'levels': set()})
            entry['n'] += 1
            if ref.get('Locator') and ref['Locator'] not in entry['locators']:
                entry['locators'].append(ref['Locator'])
            agg = (ref.get('Aggregator') or {}).get('name')
            if agg and agg not in entry['aggregators']:
                entry['aggregators'].append(agg)
            if ref.get('Provenance_Level'):
                entry['levels'].add(ref['Provenance_Level'])

        if papers:
            st.caption(f"{len(papers)} measuring paper(s) behind "
                       f"{len(selected_ligand.experimental_conditions)} condition(s)")
            for entry in sorted(papers.values(), key=lambda e: -e['n']):
                if entry['doi']:
                    st.markdown(f"📄 [{entry['doi']}]({doi_url(entry['doi'])}) "
                                f"— {entry['n']} condition(s)")
                else:
                    st.markdown(f"🔗 [{entry['url']}]({entry['url']}) "
                                f"— {entry['n']} condition(s), no DOI")
                details = []
                if entry['locators']:
                    details.append("data in " + ", ".join(entry['locators']))
                if entry['aggregators']:
                    details.append("via " + ", ".join(entry['aggregators']))
                elif 'primary' in entry['levels']:
                    details.append("transcribed directly from the paper")
                if details:
                    st.caption("; ".join(details))
        else:
            st.info("No reference information available for this ligand")


def _render_ligand_summary(selected_ligand):
    """Right detail column: selectivity trend plot, statistics, data quality."""
    st.subheader("Selectivity Trend")

    # Condition selector if multiple conditions
    condition_id = None
    if len(selected_ligand.experimental_conditions) > 1:
        st.caption(f"**{len(selected_ligand.experimental_conditions)} experimental conditions available**")

        condition_option = st.radio(
            "Display",
            ["Aggregated (median)", "Specific condition"],
            horizontal=True,
            help="Show aggregated trend across all conditions or select a specific condition"
        )

        if condition_option == "Specific condition":
            condition_ids = [c['condition_id'] for c in selected_ligand.experimental_conditions]
            condition_id = st.selectbox(
                "Select experimental condition",
                options=condition_ids,
                format_func=lambda x: x[:70] + '...' if len(x) > 70 else x
            )
    else:
        if selected_ligand.experimental_conditions:
            condition_id = selected_ligand.experimental_conditions[0]['condition_id']

    # Plot
    fig = plot_selectivity_trend(selected_ligand, condition_id)
    if fig:
        st.plotly_chart(fig, use_container_width=True)
    else:
        st.warning("Insufficient data for selectivity trend plot (need ≥2 lanthanides)")

    # Metrics
    col_m1, col_m2, col_m3, col_m4 = st.columns(4)
    col_m1.metric("Trend", selected_ligand.selectivity_trend)
    col_m2.metric("Confidence", selected_ligand.trend_confidence or 'N/A')
    sel_block = (selected_ligand.selectivity or {}) if hasattr(selected_ligand, 'selectivity') else {}
    if sel_block.get('r_squared') is not None:
        col_m3.metric("R²", f"{sel_block['r_squared']:.3f}")
    if sel_block.get('peak_element'):
        col_m3.caption(f"Peak/trough: **{sel_block['peak_element']}**")
    ci_lo = sel_block.get('slope_CI_lower')
    ci_hi = sel_block.get('slope_CI_upper')
    col_m4.metric("Slope 95% CI",
                  f"[{ci_lo:.3f}, {ci_hi:.3f}]" if (ci_lo is not None and ci_hi is not None) else "N/A")

    # Quality badge
    dq = (selected_ligand.data_quality or {}) if hasattr(selected_ligand, 'data_quality') else {}
    if dq:
        q_score = dq.get('overall_quality_score', '?')
        q_grade = dq.get('quality_grade', '?')
        n_outliers = len(dq.get('outlier_flags', []))
        badge_color = {'High': '🟢', 'Medium': '🟡', 'Low': '🔴'}.get(q_grade, '⚪')
        st.caption(
            f"{badge_color} **Data Quality**: {q_grade} (score {q_score}/100)"
            + (f" · {n_outliers} outlier flag(s) detected" if n_outliers else "")
        )


def _render_condition_explorer(selected_ligand):
    """Filterable table of every experimental condition for one ligand.

    Only rendered when the ligand has more than one condition to compare.
    """
    if len(selected_ligand.experimental_conditions) > 1:
        st.divider()
        st.subheader(f"Experimental Conditions ({len(selected_ligand.experimental_conditions)} total)")

        # Condition filters - First row: Ln conc, Acid type, Acid conc
        col_f1, col_f2, col_f3 = st.columns(3)

        with col_f1:
            # Check if Ln concentration field exists
            ln_conc_field = None
            for field_name in ['Ln_Concentration_M', 'Metal_Concentration_mM', 'Con.', 'Ln c.c./mM']:
                if selected_ligand.experimental_conditions and field_name in selected_ligand.experimental_conditions[0]:
                    ln_conc_field = field_name
                    break

            if ln_conc_field:
                ln_concs = list(set([c.get(ln_conc_field) for c in selected_ligand.experimental_conditions if c.get(ln_conc_field)]))
                ln_conc_filter = st.selectbox("Ln Concentration", ["All"] + sorted([str(x) for x in ln_concs if x]))
            else:
                st.caption("⚠️ Ln conc. not available")
                ln_conc_filter = "All"

        with col_f2:
            acids = list(set([c.get('Acid_Name', 'Unknown') for c in selected_ligand.experimental_conditions]))
            acid_filter = st.selectbox("Acid Type", ["All"] + sorted([a for a in acids if a]))

        with col_f3:
            acid_concs = list(set([c.get('Acid_Concentration_M') for c in selected_ligand.experimental_conditions if c.get('Acid_Concentration_M')]))
            if acid_concs:
                acid_conc_options = sorted([f"{x:.2f}" for x in acid_concs])
                acid_conc_filter = st.selectbox("Acid Conc (M)", ["All"] + acid_conc_options)
            else:
                st.caption("⚠️ Acid conc. not available")
                acid_conc_filter = "All"

        # Second row: Organic phase, Organic conc, Temperature
        col_f4, col_f5, col_f6 = st.columns(3)

        with col_f4:
            solvents = list(set([c.get('Solvent_Name', 'Unknown')[:30] for c in selected_ligand.experimental_conditions]))
            solvent_filter = st.selectbox("Organic Phase", ["All"] + sorted([s for s in solvents if s]))

        with col_f5:
            org_concs = list(set([c.get('Extractant_Concentration_M') for c in selected_ligand.experimental_conditions if c.get('Extractant_Concentration_M')]))
            if org_concs:
                org_conc_options = sorted([f"{x:.3f}" for x in org_concs])
                org_conc_filter = st.selectbox("Organic Conc (M)", ["All"] + org_conc_options)
            else:
                st.caption("⚠️ Org. conc. not available")
                org_conc_filter = "All"

        with col_f6:
            temps = list(set([c.get('Temperature_C') for c in selected_ligand.experimental_conditions if c.get('Temperature_C')]))
            if temps:
                temp_options = sorted([f"{x:.1f}" for x in temps])
                temp_filter = st.selectbox("Temperature (°C)", ["All"] + temp_options)
            else:
                st.caption("⚠️ Temp. not available")
                temp_filter = "All"

        # Sorting and data source filter
        st.write("")  # Spacing
        col_sort, col_source = st.columns(2)

        with col_sort:
            sort_by = st.selectbox("Sort by", ["Condition ID", "Acid concentration", "Temperature", "Elements tested"])

        # The public release carries no source labels: no filter, no column
        _meta = getattr(selected_ligand, 'metadata', None) or {}
        has_sources = bool(_meta.get('data_source') or _meta.get('data_sources'))
        source_filter = "All"
        with col_source:
            if has_sources:
                if isinstance(_meta.get('data_sources'), list):
                    source_opts = _meta['data_sources']
                else:
                    source_opts = list(set([c.get('data_source', 'Unknown') for c in selected_ligand.experimental_conditions]))
                source_filter = st.selectbox("Data Source", ["All"] + sorted(source_opts))

        # Filter conditions based on all 6 filters
        filtered_conditions = selected_ligand.experimental_conditions

        if ln_conc_filter != "All" and ln_conc_field:
            filtered_conditions = [c for c in filtered_conditions if str(c.get(ln_conc_field)) == ln_conc_filter]

        if acid_filter != "All":
            filtered_conditions = [c for c in filtered_conditions if c.get('Acid_Name') == acid_filter]

        if acid_conc_filter != "All":
            filtered_conditions = [c for c in filtered_conditions
                                 if f"{c.get('Acid_Concentration_M', 0):.2f}" == acid_conc_filter]

        if solvent_filter != "All":
            filtered_conditions = [c for c in filtered_conditions if c.get('Solvent_Name', '')[:30] == solvent_filter]

        if org_conc_filter != "All":
            filtered_conditions = [c for c in filtered_conditions
                                 if f"{c.get('Extractant_Concentration_M', 0):.3f}" == org_conc_filter]

        if temp_filter != "All":
            filtered_conditions = [c for c in filtered_conditions
                                 if f"{c.get('Temperature_C', 0):.1f}" == temp_filter]

        # Determine ligand-level data source(s) for fallback
        ligand_data_sources = []
        if hasattr(selected_ligand, 'metadata') and isinstance(selected_ligand.metadata, dict):
            if 'data_sources' in selected_ligand.metadata:
                ligand_data_sources = selected_ligand.metadata['data_sources']
            elif 'data_source' in selected_ligand.metadata:
                ligand_data_sources = [selected_ligand.metadata['data_source']]

        # For single-source ligands, use that source for all conditions
        # For multi-source ligands, conditions don't have data_source set, so filter won't work
        if source_filter != "All":
            # Only filter if ligand has single source matching filter
            if len(ligand_data_sources) == 1 and ligand_data_sources[0] == source_filter:
                pass  # Keep all conditions
            elif source_filter in ligand_data_sources:
                # Multi-source ligand - can't filter by condition without data_source tags
                pass  # Keep all conditions (filtering not possible without condition-level tags)
            else:
                filtered_conditions = []  # No conditions match

        # Create condition DataFrame with R² calculation
        cond_data = []
        for c in filtered_conditions:
            r2, slope = calculate_condition_r2(c.get('distribution_ratios', {}))

            # Format values for display (avoid Arrow serialization warnings)
            acid_conc = c.get('Acid_Concentration_M')
            acid_conc_str = f"{acid_conc:.2f}" if isinstance(acid_conc, (int, float)) else 'N/A'

            temp = c.get('Temperature_C')
            temp_str = f"{temp:.1f}" if isinstance(temp, (int, float)) else 'N/A'

            # Determine data source: use condition-level if available, else use ligand-level
            condition_source = c.get('data_source')
            if not condition_source:
                # Fallback to ligand-level metadata
                if len(ligand_data_sources) == 1:
                    condition_source = ligand_data_sources[0]
                elif len(ligand_data_sources) > 1:
                    condition_source = 'Multiple'  # Multi-source ligand
                else:
                    condition_source = 'Unknown'

            # Get DOI from metadata.references based on data source
            doi = 'N/A'
            if hasattr(selected_ligand, 'metadata') and 'references' in selected_ligand.metadata:
                for ref in selected_ligand.metadata['references']:
                    if ref.get('source_database') == condition_source:
                        doi = ref.get('doi', 'N/A')
                        break
            # Fallback to condition reference if no match found
            if doi == 'N/A':
                doi = c.get('reference', {}).get('DOI', 'N/A')
            doi_display = doi[:30] + '...' if len(doi) > 30 else doi

            # Peak element annotation (ligand-level, same for all conditions)
            _sel_block = (selected_ligand.selectivity or {}) if hasattr(selected_ligand, 'selectivity') else {}
            _peak_label = _sel_block.get('peak_element', '')

            cond_data.append({
                'Condition ID': c['condition_id'][:50] + '...' if len(c['condition_id']) > 50 else c['condition_id'],
                'Acid': c.get('Acid_Name', 'N/A'),
                'Acid Conc (M)': acid_conc_str,
                'Solvent': c.get('Solvent_Name', 'N/A')[:25],
                'Temp (°C)': temp_str,
                'Elements': len(c.get('distribution_ratios', {})),
                'Bounds (D)': condition_bounds_text(c.get('distribution_ratios', {})),
                'Co-extractants': co_extractants_text(c),
                'R²': f"{r2:.3f}" if r2 is not None else 'N/A',
                'Slope': f"{slope:.4f}" if slope is not None else 'N/A',
                'Peak': _peak_label,
                'Source': condition_source,
                'DOI': doi_display,
                'Locator': (c.get('reference') or {}).get('Locator') or ''
            })
            if not has_sources:
                del cond_data[-1]['Source']

        cond_df = pd.DataFrame(cond_data)

        # Sort
        if sort_by == "Acid concentration":
            # Sort by acid concentration, handling 'N/A'
            cond_df = cond_df.sort_values('Acid Conc (M)', key=lambda x: pd.to_numeric(x, errors='coerce').fillna(-1), ascending=False)
        elif sort_by == "Temperature":
            cond_df = cond_df.sort_values('Temp (°C)', key=lambda x: pd.to_numeric(x, errors='coerce').fillna(-1), ascending=False)
        elif sort_by == "Elements tested":
            cond_df = cond_df.sort_values('Elements', ascending=False)

        st.dataframe(cond_df, use_container_width=True, hide_index=True)

        st.caption(f"Showing {len(filtered_conditions)} of {len(selected_ligand.experimental_conditions)} conditions")

        # Add reference expander for full details
        with st.expander("📚 View Full References"):
            if filtered_conditions:
                for idx, cond in enumerate(filtered_conditions[:10]):  # Show first 10
                    ref = cond.get('reference', {})
                    if ref:
                        st.write(f"**Condition {idx+1}**: {cond['condition_id'][:60]}")
                        if ref.get('DOI'):
                            st.markdown(f"- DOI: [{ref['DOI']}]({doi_url(ref['DOI'])})")
                        elif ref.get('URL'):
                            st.markdown(f"- Source: [{ref['URL']}]({ref['URL']}) (no DOI)")
                        if ref.get('Locator'):
                            st.write(f"- Location in paper: {ref['Locator']}")
                        if ref.get('Aggregator'):
                            agg = ref['Aggregator']
                            st.write(f"- Via compilation: {agg.get('name')}"
                                     + (f" (entry {agg['entry']})" if agg.get('entry') else ""))
                        if ref.get('Notes'):
                            st.write(f"- Notes: {ref['Notes']}")
                        if ref.get('Source_Type'):
                            st.write(f"- Type: {ref['Source_Type']}")
                        if ref.get('Data_Used'):
                            st.write(f"- Data: {ref['Data_Used']}")
                        st.divider()
                if len(filtered_conditions) > 10:
                    st.caption(f"Showing 10 of {len(filtered_conditions)} references")
            else:
                st.info("No references to display")


def _render_multisource_section(selected_ligand):
    """Per-source breakdown, shown only for ligands merged from several sources."""
    # Check if ligand has multiple data sources
    is_merged = (hasattr(selected_ligand, 'metadata') and
                isinstance(selected_ligand.metadata, dict) and
                isinstance(selected_ligand.metadata.get('data_sources', []), list) and
                len(selected_ligand.metadata.get('data_sources', [])) > 1)

    if is_merged:
        st.divider()
        st.subheader("Multi-Source Comparison")
        data_sources = selected_ligand.metadata['data_sources']

        if hasattr(selected_ligand, 'selectivity') and 'trends_by_source' in selected_ligand.selectivity:
            st.write("**Selectivity Trends by Source:**")

            cols = st.columns(len(data_sources))

            for idx, (source, trend_data) in enumerate(selected_ligand.selectivity['trends_by_source'].items()):
                with cols[idx]:
                    st.metric(source, trend_data.get('trend_classification', 'N/A'))
                    if trend_data.get('r_squared'):
                        st.caption(f"R² = {trend_data['r_squared']:.3f}")
                    if trend_data.get('slope'):
                        st.caption(f"Slope = {trend_data['slope']:.4f}")

        # Show condition breakdown
        st.write("**Conditions by Source:**")
        source_counts = {}
        for c in selected_ligand.experimental_conditions:
            src = c.get('data_source', 'Unknown')
            source_counts[src] = source_counts.get(src, 0) + 1

        for src, count in sorted(source_counts.items()):
            st.write(f"- **{src}**: {count} conditions")


def _render_ligand_detail(ligands, event):
    """Render the detail view for the ligand selected in the browser table."""
    if not (event.selection and event.selection.rows):
        st.info("👆 Select a ligand from the table above to view details")
        return

    selected_idx = event.selection.rows[0]
    selected_ligand = ligands[selected_idx]

    st.divider()

    col1, col2 = st.columns([1, 2])
    with col1:
        _render_ligand_identity(selected_ligand)
    with col2:
        _render_ligand_summary(selected_ligand)

    _render_condition_explorer(selected_ligand)
    _render_multisource_section(selected_ligand)



def _render_ligand_browser_tab(db, search_text, extractant_types_filter,
                               mechanisms, trends, show_liu_only, sources,
                               required_elements, logd_min, logd_max):
    """Render the Ligand Browser: filtering, results table, and detail view."""
    ligands = _filter_ligands(db, search_text, extractant_types_filter,
                              mechanisms, trends, show_liu_only, sources,
                              required_elements, logd_min, logd_max)

    if ligands:
        event = _render_ligand_table(ligands)
        _render_ligand_detail(ligands, event)
    else:
        st.warning("No ligands found matching your filter criteria. Try adjusting the filters in the sidebar.")



def main():
    # Load database
    try:
        db = load_database()
    except Exception as e:
        st.error(f"Error loading database: {e}")
        st.info("Make sure `data/sx_ligands_merged_v14_flat.jsonl` exists.")
        return

    # Header
    st.title("🧪 REEL-QUEST")
    st.markdown("**Rare-Earth Extraction Ligands - Quantitative Understanding of Elemental Selectivity Trends**")
    st.markdown("Interactive dashboard for browsing and visualizing ligand extraction data with uncertainty quantification")
    if get_config("RELQUEST_DB"):
        st.warning(f"Previewing `{db.db_path}` (RELQUEST_DB), not the canonical database.")

    (page, show_liu_only, search_text, extractant_types_filter, mechanisms,
     trends, sources, required_elements, logd_min, logd_max) = _render_sidebar(db)
    #=========================================================================
    # AI Assistant Page (early return — bypasses database content)
    #=========================================================================
    if page == "🤖 AI Assistant":
        _render_ai_assistant_tab(db)
        return
    all_ligands = _render_overview_tab(db, show_liu_only)

    _render_single_element_tab(db, all_ligands)

    _render_logsf_tab(db, all_ligands)

    _render_ligand_browser_tab(db, search_text, extractant_types_filter,
                               mechanisms, trends, show_liu_only, sources,
                               required_elements, logd_min, logd_max)



if __name__ == "__main__":
    main()
