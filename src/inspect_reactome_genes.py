
from collections import Counter

import requests

PATHWAY_ID = "R-HSA-913531"
url = (
    "https://reactome.org/ContentService/data/"
    f"participants/{PATHWAY_ID}"
)

response = requests.get(url, timeout=60)
response.raise_for_status()
records = response.json()

print(f"Participant records: {len(records)}")

print("\n--- SPECIES ---")
print(Counter(
    record.get("speciesName", "Missing")
    for record in records
))

print("\n--- ENTITY CLASSES ---")
print(Counter(
    record.get("schemaClass", "Missing")
    for record in records
))

print("\n--- EXAMPLE RECORDS ---")
for record in records[:5]:
    print({
        "displayName": record.get("displayName"),
        "speciesName": record.get("speciesName"),
        "schemaClass": record.get("schemaClass"),
        "refEntities": record.get("refEntities"),
    })
