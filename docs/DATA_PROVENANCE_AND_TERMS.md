# Data provenance and third-party terms

Where the values in this database come from, on what basis they are
redistributed, and how to credit each source.

Every condition names the paper that measured its values (`reference.DOI`) and
the figure, table or text page it was read from (`reference.Locator`), and
every condition was checked against that paper. What this project
redistributes is transcribed numerical results with that attribution, plus its
own curation and derived fields (see `LICENSE-DATA`).

## Sources

The values were compiled from the following sources and from papers read
directly for this project.

### Liu et al. 2022, JACS Au

- **Citation:** Liu, T.; Johnson, K. R.; Jansone-Popova, S.; Jiang, D.
  "Advancing Rare-Earth Separation by Machine Learning." *JACS Au* **2022**, 2,
  1428–1434. https://doi.org/10.1021/jacsau.2c00122
- **Taken:** distribution ratios from the Supporting Information, a
  compilation of 45 source papers; each condition cites its source paper.
- **Terms basis:** the article is open access under CC BY; the Supporting
  Information is published with it.

### SAFE — Separation Archive for Elements

- **Citation (requested by SAFE):** SAFE asks users to "cite the following
  paper if you use data from this website": Garcia Leite, L. da; Zhang, B.;
  Acar, Z.; Elowitt, J.; Augustine, L. J.; Clark, A. E.; Karamalis, V.;
  Perez, D.; Schrier, J.; Taylor, M.; Yang, P. "Creation of the Separation
  Archive for Elements (SAFE) Database." *Solvent Extraction and Ion Exchange*
  **2025**, 43(7), 847–852. https://doi.org/10.1080/07366299.2025.2564381
- **Taken:** lanthanide, Y and Sc distribution ratios from a SAFE export
  (2026-01-26) for extractants with a valid structure, split into one
  condition per experiment and checked against the source papers. The raw
  export is not redistributed.
- **Terms basis:** SAFE describes itself as "a centralized, publicly accessible
  database" of "AI-Ready, Machine Readable" data and asks for the citation
  above. Only transcribed, attributed values are redistributed.

### OECD-NEA IDEaL — International Database for Extractant Ligands

- **Citation:** OECD Nuclear Energy Agency, IDEaL — International Database for
  Extractant Ligands, https://www.oecd-nea.org/ideal/
- **Taken:** distribution ratios with the reference each IDEaL entry cites,
  checked against those papers.

### Papers read directly

Values transcribed by this project from the paper's tables or from
calibrated figures, including Pramanik, S. *et al.* *J. Am. Chem. Soc.*
**2024**, 146, 25669–25679 (https://doi.org/10.1021/jacs.4c07332, Supporting
Information Tables S1 and S2). Each is cited in its conditions' `reference.DOI`.

## How to credit

- The database: `CITATION.cff` (GitHub's "Cite this repository").
- Specific measurements: the papers in their `reference.DOI`, all of them in
  `paper/references.bib` (BibTeX).
- Liu et al. 2022 and IDEaL: as above, when your work relies on the compiled
  data as a whole.
