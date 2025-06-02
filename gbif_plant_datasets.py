import requests
import random

# Constants for GBIF API
GBIF_DATASET_BASE = "https://api.gbif.org/v1/dataset/search"
GBIF_OCCURRENCE_BASE = "https://api.gbif.org/v1/occurrence/search"
PLANT_TAXON_KEY = 6  # Key for Plantae in GBIF backbone taxonomy
COMMERCIAL_LICENSES = ["CC0_1_0", "CC_BY_4_0"]

# List of geospatial and taxonomic issues from GBIF documentation
GEOSPATIAL_ISSUES = [
    "ZERO_COORDINATE", "COORDINATE_OUT_OF_RANGE", "COUNTRY_COORDINATE_MISMATCH", "COORDINATE_INVALID",
    "GEODETIC_DATUM_ASSUMED_WGS84", "GEODETIC_DATUM_INVALID", "FOOTPRINT_SRS_INVALID", "FOOTPRINT_WKT_MISMATCH",
    "FOOTPRINT_WKT_INVALID", "COUNTRY_MISMATCH", "COUNTRY_DERIVED_FROM_COORDINATES", "COUNTRY_INVALID",
    "CONTINENT_COORDINATE_MISMATCH", "CONTINENT_DERIVED_FROM_COUNTRY", "CONTINENT_INVALID", "COORDINATE_ROUNDED",
    "COORDINATE_REPROJECTED", "COORDINATE_REPROJECTION_SUSPICIOUS", "COORDINATE_REPROJECTION_FAILED",
    "COORDINATE_UNCERTAINTY_METERS_INVALID"
]
TAXONOMIC_ISSUES = [
    "TAXON_MATCH_NONE", "TAXON_MATCH_HIGHERRANK", "TAXON_MATCH_FUZZY", "TAXON_MATCH_AGGREGATE",
    "TAXON_MATCH_PARENT", "TAXON_MATCH_ACCEPTED", "TAXON_MATCH_SYNONYM", "TAXON_MATCH_AMBIGUOUS"
]

def get_first_page_dataset_keys(license_types, limit=200):
    keys = []
    for license_type in license_types:
        params = {
            "type": "OCCURRENCE",
            "license": license_type,
            "limit": limit,
            "offset": 0
        }
        resp = requests.get(GBIF_DATASET_BASE, params=params)
        resp.raise_for_status()
        data = resp.json()
        for ds in data.get("results", []):
            keys.append(ds["key"])
    return keys

def has_plantae_occurrence(dataset_key):
    params = {
        "datasetKey": dataset_key,
        "taxonKey": PLANT_TAXON_KEY,
        "limit": 1
    }
    resp = requests.get(GBIF_OCCURRENCE_BASE, params=params)
    resp.raise_for_status()
    return resp.json().get("count", 0) > 0

def analyze_occurrence_quality(occurrences):
    total = 0
    fit = 0
    for idx, rec in enumerate(occurrences, 1):
        if rec.get("kingdomKey") != PLANT_TAXON_KEY:
            continue
        total += 1
        issues = rec.get("issues", [])
        if not any(issue in GEOSPATIAL_ISSUES or issue in TAXONOMIC_ISSUES for issue in issues):
            fit += 1
        if total > 0 and total % 50 == 0:
            print(f"Record {total}: Issues: {issues}")
    percent = (fit / total) * 100 if total > 0 else 0.0
    print(f"\nOut of {total} plant records (kingdomKey==6), {fit} have no geospatial or taxonomic issues.")
    print(f"Percentage of fit-to-go records: {percent:.2f}%")

def get_plant_occurrences_for_dataset(dataset_key, n=1000, page_limit=300):
    occurrences = []
    offset = 0
    while len(occurrences) < n:
        params = {
            "datasetKey": dataset_key,
            "taxonKey": PLANT_TAXON_KEY,
            "hasCoordinate": "true",
            "limit": min(page_limit, n - len(occurrences)),
            "offset": offset
        }
        resp = requests.get(GBIF_OCCURRENCE_BASE, params=params)
        resp.raise_for_status()
        data = resp.json()
        records = data.get("results", [])
        if not records:
            break
        occurrences.extend(records)
        if len(records) < page_limit:
            break
        offset += page_limit
    return occurrences[:n]

def get_occurrence_count(dataset_key, taxon_key=None):
    params = {
        "datasetKey": dataset_key,
        "limit": 0
    }
    if taxon_key is not None:
        params["taxonKey"] = taxon_key
    resp = requests.get(GBIF_OCCURRENCE_BASE, params=params)
    resp.raise_for_status()
    return resp.json().get("count", 0)

if __name__ == "__main__":
    print("Filtering first page of datasets at the dataset level (license/type)...")
    dataset_keys = get_first_page_dataset_keys(COMMERCIAL_LICENSES, limit=200)
    print(f"Number of occurrence datasets in first page: {len(dataset_keys)}")
    if not dataset_keys:
        print("No datasets found.")
    else:
        print("\nCalculating Plantae percentage for each dataset...")
        plantae_percentages = {}
        for key in dataset_keys:
            total_count = get_occurrence_count(key)
            if total_count == 0:
                continue
            plantae_count = get_occurrence_count(key, PLANT_TAXON_KEY)
            if plantae_count == 0:
                continue
            percentage = (plantae_count / total_count) * 100
            plantae_percentages[key] = percentage
        print(f"Number of datasets with at least one Plantae record: {len(plantae_percentages)}")
        print("Plantae percentages for all datasets with at least one Plantae record:")
        for k, v in plantae_percentages.items():
            print(f"{k}: {v:.2f}%")
        # Only keep datasets with >50% Plantae
        high_plantae_datasets = [k for k, v in plantae_percentages.items() if v > 50]
        print(f"\nNumber of datasets with >50% Plantae: {len(high_plantae_datasets)}")
        if high_plantae_datasets:
            sample = random.sample(high_plantae_datasets, min(10, len(high_plantae_datasets)))
            for i, selected_dataset in enumerate(sample, 1):
                print(f"\nAnalyzing quality for dataset {i}: {selected_dataset}")
                plant_occurrences = get_plant_occurrences_for_dataset(selected_dataset, n=1000)
                if not plant_occurrences:
                    print("No Plantae records with coordinates found in this dataset.")
                    continue
                num_kingdom_plants = sum(1 for rec in plant_occurrences if rec.get("kingdomKey") == PLANT_TAXON_KEY)
                print(f"Number of records fetched: {len(plant_occurrences)}; Number with kingdomKey==6: {num_kingdom_plants}")
                analyze_occurrence_quality(plant_occurrences)
        else:
            print("No datasets with >50% Plantae found.") 