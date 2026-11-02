"""Which repository this checkout is.

'dev'     the developer repository: full provenance (data sources, compilation
          references, curation history) and the pipelines that built the data.
'public'  the published repository, exported by scripts/export_public_repo.py,
          which replaces this file. Records there carry only
          rel_quest.schema.PUBLIC_FIELDS.
"""

PROFILE = "public"
