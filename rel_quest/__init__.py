"""
REEL-QUEST - Rare-Earth Extraction Ligands
(Quantitative Understanding of Elemental Selectivity Trends)
A Python library for querying rare earth ligand extraction data
"""

__version__ = '1.0.0'
__description__ = 'Rare earth ligand database with mechanism classification and selectivity trends'

# Main exports
from .schema import SCHEMA_VERSION, KNOWN_SCHEMA_VERSIONS
from .core import DatabaseLoadError, LigandDatabase, LigandQuery
from .utils import load_vocabulary, canonicalize_smiles, validate_smiles

__all__ = [
    'SCHEMA_VERSION',
    'KNOWN_SCHEMA_VERSIONS',
    'LigandDatabase',
    'LigandQuery',
    'DatabaseLoadError',
    'load_vocabulary',
    'canonicalize_smiles',
    'validate_smiles'
]
