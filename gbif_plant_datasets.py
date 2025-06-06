import requests
import random
import matplotlib.pyplot as plt
import numpy as np
import os
from collections import defaultdict, Counter
from matplotlib.ticker import FuncFormatter

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

def get_dataset_keys(license_types, page_size=500, datasets_per_license=500):
    """
    Get dataset keys from the GBIF API by selecting one random page of datasets for each license type.
    
    Args:
        license_types: List of license types to include
        page_size: Number of results per page request (max 1000 for GBIF API)
        datasets_per_license: Target number of datasets to retrieve per license type
        
    Returns:
        List of dataset keys
    """
    keys = []
    
    for license_type in license_types:
        # Get a count of total datasets for this license type
        count_params = {
            "type": "OCCURRENCE",
            "license": license_type,
            "limit": 0  # Just get the count, not actual results
        }
        try:
            count_resp = requests.get(GBIF_DATASET_BASE, params=count_params)
            count_resp.raise_for_status()
            count_data = count_resp.json()
            total_count = count_data.get("count", 0)
            
            if total_count == 0:
                print(f"No datasets found for license type: {license_type}")
                continue
                
            # Calculate how many pages we can sample from with page_size
            max_page = (total_count // page_size) + 1
            print(f"Found {total_count} datasets for license type {license_type} across {max_page} pages")
            
            # If there are fewer datasets than our page_size, adjust to get all available
            if total_count <= page_size:
                print(f"Total datasets for {license_type} is less than page size. Retrieving all available datasets.")
                random_page = 1
                actual_page_size = total_count
            else:
                # Randomly select one page
                random_page = random.randint(1, max_page)
                actual_page_size = min(page_size, datasets_per_license)
                
            # Calculate offset based on the selected page
            offset = (random_page - 1) * page_size
                
            params = {
                "type": "OCCURRENCE",
                "license": license_type,
                "limit": actual_page_size,
                "offset": offset
            }
            
            print(f"Requesting page {random_page} with {actual_page_size} datasets for license {license_type}")
            
            resp = requests.get(GBIF_DATASET_BASE, params=params)
            resp.raise_for_status()
            data = resp.json()
            results = data.get("results", [])
            
            # Add the new keys
            new_keys = [ds["key"] for ds in results]
            keys.extend(new_keys)
            
            print(f"Retrieved {len(new_keys)} datasets for license {license_type}")
                
        except Exception as e:
            print(f"Error retrieving datasets for license type {license_type}: {e}")
    
    # Remove any duplicates that might have occurred
    unique_keys = list(set(keys))
    print(f"Retrieved {len(unique_keys)} unique datasets out of {len(keys)} total fetched")
    
    return unique_keys

def count_issues_in_records(occurrences, verbose=True):
    """
    Count issues in a collection of occurrence records and return detailed statistics.
    
    Args:
        occurrences: List of occurrence records
        verbose: Whether to print progress and summary
        
    Returns:
        Dictionary containing record counts and issue statistics
    """
    total = 0
    fit = 0
    issue_counts = {issue: 0 for issue in GEOSPATIAL_ISSUES + TAXONOMIC_ISSUES}
    
    for idx, rec in enumerate(occurrences, 1):
        if rec.get("kingdomKey") != PLANT_TAXON_KEY:
            continue
            
        total += 1
        issues = rec.get("issues", [])
        
        # Count each issue type
        for issue in issues:
            if issue in issue_counts:
                issue_counts[issue] += 1
        
        # Check if record has no issues
        if not any(issue in GEOSPATIAL_ISSUES or issue in TAXONOMIC_ISSUES for issue in issues):
            fit += 1
            
        # Print progress if verbose
        if verbose and total > 0 and total % 50 == 0:
            print(f"Record {total}: Issues: {issues}")
    
    # Calculate percentage of fit records
    percent = (fit / total) * 100 if total > 0 else 0.0
    
    # Only return issues that actually occurred
    active_issues = {k: v for k, v in issue_counts.items() if v > 0}
    
    # Print summary if verbose
    if verbose and total > 0:
        print(f"\nOut of {total} plant records (kingdomKey==6), {fit} have no geospatial or taxonomic issues.")
        print(f"Percentage of fit-to-go records: {percent:.2f}%")
    
    return {
        "total_records": total,
        "fit_records": fit,
        "fit_percent": percent,
        "issue_counts": active_issues
    }

def analyze_occurrence_quality(occurrences):
    """
    Analyze quality of occurrence records and print summary statistics.
    
    Args:
        occurrences: List of occurrence records
        
    Returns:
        Tuple of (total_records, fit_records, fit_percent)
    """
    results = count_issues_in_records(occurrences)
    return results["total_records"], results["fit_records"], results["fit_percent"]

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
    """
    Get count of occurrences for a dataset, optionally filtered by taxon key
    
    Args:
        dataset_key: GBIF dataset key
        taxon_key: Optional taxonomy key to filter by (e.g., PLANT_TAXON_KEY)
        
    Returns:
        Count of matching occurrences
    """
    params = {
        "datasetKey": dataset_key,
        "limit": 0
    }
    if taxon_key is not None:
        params["taxonKey"] = taxon_key
    resp = requests.get(GBIF_OCCURRENCE_BASE, params=params)
    resp.raise_for_status()
    return resp.json().get("count", 0)

def calculate_plantae_percentages(dataset_keys):
    """
    Calculate percentage of Plantae records for each dataset
    
    Args:
        dataset_keys: List of GBIF dataset keys
        
    Returns:
        Dictionary mapping dataset keys to their Plantae percentages
    """
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
        
    return plantae_percentages

def create_plantae_percentage_cdf(plantae_percentages):
    """
    Create a cumulative distribution function (CDF) graph showing the number of datasets
    that have a plantae percentage greater than or equal to each threshold.
    
    Args:
        plantae_percentages: Dictionary mapping dataset keys to their Plantae percentages
        
    Returns:
        None (saves CDF graph to file)
    """
    if not plantae_percentages:
        print("No plantae percentages to plot.")
        return
        
    # Get all percentage values
    percentages = list(plantae_percentages.values())
    
    # Create percentage thresholds from 0 to 100
    thresholds = list(range(0, 101, 5))  # 0, 5, 10, ..., 100
    
    # Count datasets above each threshold
    counts = [sum(1 for pct in percentages if pct >= threshold) for threshold in thresholds]
    
    # Create the figure with better proportions
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Plot the CDF with better styling
    ax.plot(thresholds, counts, marker='o', markersize=6, linestyle='-', 
            linewidth=2.5, color='forestgreen', zorder=3)
    
    # Fill the area under the curve for better visualization
    ax.fill_between(thresholds, counts, alpha=0.2, color='forestgreen')
    
    # Add a vertical line at 50% threshold
    ax.axvline(x=50, color='red', linestyle='--', alpha=0.7, 
              linewidth=2, label='Current threshold (50%)', zorder=2)
    
    # Add data point at 50% for clarity
    count_at_50 = sum(1 for pct in percentages if pct >= 50)
    ax.plot(50, count_at_50, 'ro', ms=10, zorder=4)
    
    # Annotate the point at 50% with a clearer label
    ax.annotate(f'{count_at_50} datasets',
                xy=(50, count_at_50),
                xytext=(58, count_at_50 + 2),  # Position label away from the line
                arrowprops=dict(facecolor='black', shrink=0.05, width=1.5, 
                               headwidth=8, alpha=0.7),
                fontsize=11,
                fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8))
    
    # Customize the plot with better formatting
    ax.set_xlabel('Plantae Percentage Threshold (%)', fontsize=12)
    ax.set_ylabel('Number of Datasets with ≥ X% Plantae Records', fontsize=12)
    ax.set_title('Cumulative Distribution of Datasets by Plantae Percentage', 
                fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, linestyle='--', alpha=0.7, zorder=1)
    
    # Add x-axis ticks at every 10%
    ax.set_xticks(range(0, 101, 10))
    
    # Set y-axis to start from 0 and set a reasonable top limit
    total_datasets = len(percentages)
    ax.set_ylim(bottom=0, top=total_datasets * 1.1)
    
    # Show the total number of datasets at 0% threshold
    ax.annotate(f'Total: {total_datasets} datasets',
                xy=(5, total_datasets),
                xytext=(5, total_datasets * 0.95),
                fontsize=10,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.7))
    
    # Add legend with better formatting
    ax.legend(loc='upper right', frameon=True, framealpha=0.8)
    
    # Save the figure
    save_figure(fig, 'plantae_percentage_cdf.png')
    plt.close()

def analyze_quality_issues_across_datasets(datasets_list, n_records=1000):
    """
    Analyze quality issues across multiple datasets and create a bar graph
    showing the percentage of occurrences affected by each issue flag.
    
    Args:
        datasets_list: List of dataset keys to analyze
        n_records: Number of records to fetch per dataset
    
    Returns:
        None (saves bar graph to file)
    """
    # Combined issue counts across all datasets
    combined_issue_counts = {}
    total_records = 0
    
    print("\nAnalyzing quality issues across all datasets...")
    
    # Process each dataset
    for i, dataset_key in enumerate(datasets_list, 1):
        print(f"Processing dataset {i}/{len(datasets_list)}: {dataset_key}")
        occurrences = get_plant_occurrences_for_dataset(dataset_key, n=n_records)
        if not occurrences:
            print("No Plantae records with coordinates found in this dataset.")
            continue
        
        # Use our common function to count issues
        results = count_issues_in_records(occurrences, verbose=False)
        
        # Aggregate the counts
        for issue, count in results["issue_counts"].items():
            combined_issue_counts[issue] = combined_issue_counts.get(issue, 0) + count
        
        total_records += results["total_records"]
        print(f"Processed {results['total_records']} plant records from dataset {dataset_key}")
    
    # Calculate percentages
    if total_records > 0:
        issue_percentages = {issue: (count / total_records) * 100 
                            for issue, count in combined_issue_counts.items()}
        
        # Create the bar graph
        create_issue_bar_graph(issue_percentages, total_records)
    else:
        print("No plant records found across all datasets.")

def save_figure(fig, filename, dpi=300):
    """
    Save a matplotlib figure with consistent settings
    
    Args:
        fig: The matplotlib figure to save
        filename: The filename (without path)
        dpi: Resolution for the saved image
    """
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), filename)
    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
    print(f"Figure saved to: {output_path}")
    return output_path

def create_issue_bar_graph(issue_percentages, total_records):
    """
    Create a bar graph showing the percentage of occurrences affected by each issue flag.
    
    Args:
        issue_percentages: Dictionary mapping issue flags to their percentages
        total_records: Total number of records analyzed
    
    Returns:
        None (saves bar graph to file)
    """
    # Sort issues by percentage (descending)
    sorted_issues = sorted(issue_percentages.items(), key=lambda x: x[1], reverse=True)
    
    # Extract issue names and percentages
    issues = [issue for issue, _ in sorted_issues]
    percentages = [pct for _, pct in sorted_issues]
    
    # Create figure with appropriate size based on number of issues
    fig, ax = plt.subplots(figsize=(12, max(8, len(issues) * 0.4)))
    
    # Create horizontal bar chart
    bars = ax.barh(issues, percentages, color='skyblue')
    
    # Add percentage labels to the bars
    for bar in bars:
        width = bar.get_width()
        ax.text(width + 0.5, bar.get_y() + bar.get_height()/2, 
                f'{width:.2f}%', ha='left', va='center')
    
    # Customize the plot
    ax.set_xlabel('Percentage of Records (%)')
    ax.set_title(f'Quality Issue Flags Distribution\n(Based on {total_records} plant records from GBIF)')
    ax.grid(axis='x', linestyle='--', alpha=0.7)
    
    # Adjust layout and save
    plt.tight_layout()
    save_figure(fig, 'quality_issues_bar_graph.png')
    
    # Create pie charts for geospatial vs taxonomic issues
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Separate geospatial and taxonomic issues
    geo_issues = {k: v for k, v in issue_percentages.items() if k in GEOSPATIAL_ISSUES}
    taxo_issues = {k: v for k, v in issue_percentages.items() if k in TAXONOMIC_ISSUES}
    
    # Sort issues by value for better visualization
    sorted_geo_issues = dict(sorted(geo_issues.items(), key=lambda x: x[1], reverse=True))
    sorted_taxo_issues = dict(sorted(taxo_issues.items(), key=lambda x: x[1], reverse=True))
    
    # Function to make labels more readable
    def shorten_labels(k):
        return k.replace('_', ' ').title()
    
    # Create two subplots for comparison
    if geo_issues:
        # Create pie chart with better label spacing
        wedges, texts, autotexts = ax1.pie(
            sorted_geo_issues.values(), 
            autopct='%1.1f%%',
            textprops={'fontsize': 9},
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
            pctdistance=0.85
        )
        
        # Place labels outside of pie with lines connecting to wedges
        ax1.legend(
            wedges,
            [shorten_labels(k) for k in sorted_geo_issues.keys()],
            title="Geospatial Issues",
            loc="center left",
            bbox_to_anchor=(0.9, 0, 0.5, 1)
        )
        
        ax1.set_title('Geospatial Issues', pad=20)
    else:
        ax1.text(0.5, 0.5, 'No geospatial issues found', ha='center', va='center')
    
    if taxo_issues:
        # Create pie chart with better label spacing
        wedges, texts, autotexts = ax2.pie(
            sorted_taxo_issues.values(), 
            autopct='%1.1f%%',
            textprops={'fontsize': 9},
            wedgeprops={'linewidth': 1, 'edgecolor': 'white'},
            pctdistance=0.85
        )
        
        # Place labels in a legend instead of directly on the pie
        ax2.legend(
            wedges,
            [shorten_labels(k) for k in sorted_taxo_issues.keys()],
            title="Taxonomic Issues",
            loc="center right",
            bbox_to_anchor=(0.1, 0, 0.5, 1)
        )
        
        ax2.set_title('Taxonomic Issues', pad=20)
    else:
        ax2.text(0.5, 0.5, 'No taxonomic issues found', ha='center', va='center')
    
    plt.suptitle('Distribution of Geospatial vs Taxonomic Issues')
    plt.tight_layout()
    
    # Save the comparison chart
    save_figure(fig, 'issue_types_comparison.png')
    
    # Close the plots to free memory
    plt.close('all')

def analyze_dataset_issue_profile(dataset_key, n_records=1000):
    """
    Analyze the issue profile for a specific dataset.
    
    Args:
        dataset_key: The GBIF dataset key
        n_records: Number of records to fetch
        
    Returns:
        Dictionary with issue counts and total records
    """
    # Get occurrences
    occurrences = get_plant_occurrences_for_dataset(dataset_key, n=n_records)
    if not occurrences:
        print(f"No Plantae records with coordinates found in dataset {dataset_key}.")
        return None
    
    # Use our common function to count issues
    results = count_issues_in_records(occurrences, verbose=False)
    
    return {
        "dataset_key": dataset_key,
        "total_records": results["total_records"],
        "issue_counts": results["issue_counts"]
    }

def create_dataset_comparison_chart(dataset_profiles):
    """
    Create a chart comparing issue profiles across datasets.
    
    Args:
        dataset_profiles: List of dataset issue profiles
        
    Returns:
        None (saves chart to file)
    """
    if not dataset_profiles:
        print("No valid dataset profiles to compare.")
        return
    
    # Get all unique issues across all datasets
    all_issues = set()
    for profile in dataset_profiles:
        all_issues.update(profile["issue_counts"].keys())
    
    # Sort issues by frequency across all datasets
    issue_total_counts = {issue: sum(profile["issue_counts"].get(issue, 0) 
                                    for profile in dataset_profiles) 
                         for issue in all_issues}
    sorted_issues = sorted(issue_total_counts.items(), key=lambda x: x[1], reverse=True)
    sorted_issue_names = [issue for issue, _ in sorted_issues]
    
    # Prepare data for plotting
    dataset_keys = [profile["dataset_key"] for profile in dataset_profiles]
    dataset_labels = [f"DS{i+1}" for i in range(len(dataset_keys))]
    
    # Map actual keys to labels
    key_to_label = dict(zip(dataset_keys, dataset_labels))
    label_to_key = dict(zip(dataset_labels, dataset_keys))
    
    # Print the mapping for reference
    print("\nDataset key to label mapping:")
    for key, label in key_to_label.items():
        print(f"{label}: {key}")
    
    # Calculate percentages for each dataset
    data = []
    for profile in dataset_profiles:
        total = profile["total_records"]
        if total == 0:
            continue
            
        percentages = [(profile["issue_counts"].get(issue, 0) / total) * 100 
                      for issue in sorted_issue_names]
        data.append(percentages)
    
    # Create a heatmap if we have multiple datasets
    if len(data) > 1 and len(sorted_issue_names) > 0:
        # Create figure with enough space for long issue names
        fig_width = max(16, len(sorted_issue_names) * 1.0)
        fig_height = max(10, len(dataset_labels) * 0.8)
        fig, ax = plt.subplots(figsize=(fig_width, fig_height))
        
        # Shorten issue names for better display
        shortened_issue_names = [issue.replace('_', '\n') for issue in sorted_issue_names]
        
        # Create heatmap
        im = ax.imshow(data, aspect='auto', cmap='YlOrRd')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Percentage of Records (%)')
        
        # Add labels
        ax.set_yticks(range(len(dataset_labels)))
        ax.set_yticklabels(dataset_labels)
        ax.set_xticks(range(len(shortened_issue_names)))
        ax.set_xticklabels(shortened_issue_names, rotation=45, ha='right', fontsize=9)
        
        # Set tick parameters
        ax.tick_params(axis='x', which='major', pad=5)
        
        # Create padding between labels and plot
        plt.subplots_adjust(bottom=0.2)
        
        ax.set_xlabel('Issue Type', fontsize=12, labelpad=10)
        ax.set_ylabel('Dataset', fontsize=12)
        ax.set_title('Percentage of Records Affected by Each Issue Type Across Datasets', 
                    fontsize=14, pad=20)
        
        # Add text annotations, but only if the percentage is significant
        for i in range(len(dataset_labels)):
            for j in range(len(sorted_issue_names)):
                if j < len(data[i]):
                    if data[i][j] > 1.0:  # Only show text if percentage > 1%
                        text = f"{data[i][j]:.1f}%" 
                        text_color = "black" if data[i][j] < 50 else "white"
                        ax.text(j, i, text, ha="center", va="center", 
                                color=text_color, fontsize=8)
        
        plt.tight_layout()
        save_figure(fig, 'dataset_issue_comparison.png')
        print("\nDataset comparison heatmap saved.")
    
    # Create a summary radar chart for the top 5 issues
    if sorted_issue_names and len(data) > 0:
        # Select top 5 issues (or all if less than 5)
        top_issues = sorted_issue_names[:min(5, len(sorted_issue_names))]
        
        # Make issue labels more readable for display
        readable_issues = [issue.replace('_', ' ').title() for issue in top_issues]
        
        # Create radar chart with more space
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, polar=True)
        
        # Number of variables
        N = len(top_issues)
        
        # Compute angle for each axis
        angles = [n / float(N) * 2 * np.pi for n in range(N)]
        angles += angles[:1]  # Close the loop
        
        # Use different colors for each dataset
        colors = plt.cm.tab10(np.linspace(0, 1, len(dataset_labels)))
        
        # Draw the chart for each dataset
        for i, (label, percentages) in enumerate(zip(dataset_labels, data)):
            # Select only the top issues
            values = [percentages[sorted_issue_names.index(issue)] for issue in top_issues]
            values += values[:1]  # Close the loop
            
            # Plot values with distinct colors and better styling
            ax.plot(angles, values, linewidth=2.5, linestyle='solid', label=label, color=colors[i])
            ax.fill(angles, values, alpha=0.25, color=colors[i])
        
        # Fix axis to go in the right order and start at 12 o'clock
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        
        # Draw axis lines for each angle and label with better spacing
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(readable_issues, fontsize=10)
        
        # Add gridlines with better visibility
        ax.grid(True, linestyle='--', alpha=0.7, linewidth=0.5)
        
        # Draw y-axis labels (percentages) with better spacing
        ax.set_rlabel_position(0)
        max_value = max([max(row) for row in data]) if data else 0
        # Ensure we have enough tick marks but not too many
        num_ticks = min(5, int(max_value/20)+1)
        tick_values = [20, 40, 60, 80, 100][:num_ticks] if num_ticks > 0 else [20]
        plt.yticks(tick_values, 
                   [f"{v}%" for v in tick_values], 
                   color="grey", size=9)
        
        # Set reasonable y-axis limit
        plt.ylim(0, max(100, min(max_value * 1.2, 120)))
        
        # Add legend in a better position with better formatting
        legend = plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1),
                          frameon=True, framealpha=0.8, fontsize=9)
        legend.set_title("Datasets", prop={'size': 10, 'weight': 'bold'})
        
        plt.title('Top Issues Comparison Across Datasets', size=16, y=1.1, fontweight='bold')
        
        # Save the radar chart
        save_figure(fig, 'top_issues_radar.png')
        print("Top issues radar chart saved.")
    
    plt.close('all')

def get_recently_modified_datasets(days=30, license_types=None, limit=500):
    """
    Find datasets that were modified in the last specified number of days.
    
    Args:
        days: Number of days to look back for modifications
        license_types: List of license types to filter by (e.g., COMMERCIAL_LICENSES)
        limit: Maximum number of datasets to return per license type
        
    Returns:
        List of dataset keys that were modified in the specified time period
    """
    from datetime import datetime, timedelta
    
    # Calculate the date from days ago
    today = datetime.now()
    from_date = today - timedelta(days=days)
    from_date_str = from_date.strftime("%Y-%m-%d")
    
    print(f"\nFinding datasets modified since {from_date_str}...")
    
    modified_keys = []
    total_fetched = 0
    
    # If no specific license types provided, use a default empty list for one iteration
    license_list = license_types if license_types else [None]
    
    for license_type in license_list:
        offset = 0
        more_results = True
        
        while more_results and len(modified_keys) < limit:
            params = {
                "type": "OCCURRENCE",
                "limit": 100,
                "offset": offset
            }
            
            if license_type:
                params["license"] = license_type
                
            resp = requests.get(GBIF_DATASET_BASE, params=params)
            resp.raise_for_status()
            data = resp.json()
            
            if not data.get("results"):
                more_results = False
                break
                
            # Check each dataset's modified date
            for ds in data.get("results", []):
                total_fetched += 1
                
                if "modified" in ds:
                    modified_date_str = ds["modified"]
                    # Parse the date from ISO format (handling possible timezone info)
                    modified_date = modified_date_str.split('T')[0]  # Extract just the date part
                    modified_date = datetime.strptime(modified_date, "%Y-%m-%d")
                    
                    if modified_date >= from_date:
                        modified_keys.append(ds["key"])
                        
            offset += len(data.get("results", []))
            
            # Check if we've reached the end of results
            if data.get("endOfRecords", True):
                more_results = False
                
    print(f"Found {len(modified_keys)} datasets modified in the last {days} days (checked {total_fetched} datasets)")
    return modified_keys

def analyze_recently_modified_datasets(days=30, license_types=None):
    """
    Analyze datasets that were modified recently and plot statistics about plant instances.
    
    Args:
        days: Number of days to look back
        license_types: List of license types to filter by
        
    Returns:
        None (generates and saves a graph)
    """
    # Get recently modified datasets
    recent_dataset_keys = get_recently_modified_datasets(days, license_types)
    
    if not recent_dataset_keys:
        print("No datasets were modified in the specified time period.")
        return
        
    # Calculate plantae percentages for these datasets
    plantae_percentages = calculate_plantae_percentages(recent_dataset_keys)
    
    if not plantae_percentages:
        print("None of the recently modified datasets contain Plantae records.")
        return
        
    # Categorize datasets by plantae percentage
    ranges = [(0, 25), (25, 50), (50, 75), (75, 100)]
    range_counts = {f"{r[0]}-{r[1]}%": 0 for r in ranges}
    
    for _, percentage in plantae_percentages.items():
        for lower, upper in ranges:
            if lower <= percentage < upper or (upper == 100 and percentage == 100):
                range_counts[f"{lower}-{upper}%"] += 1
                break
    
    # Create a bar chart
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Get the ranges and counts
    categories = list(range_counts.keys())
    counts = list(range_counts.values())
    
    # Create the bar chart with enhanced styling
    bars = ax.bar(categories, counts, color='lightgreen', edgecolor='darkgreen', alpha=0.8)
    
    # Add count labels on top of each bar
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                   f'{int(height)}', ha='center', va='bottom', fontweight='bold')
    
    # Calculate the total number of recently modified datasets with plant records
    total_with_plants = sum(counts)
    total_recent = len(recent_dataset_keys)
    plant_percentage = (total_with_plants / total_recent * 100) if total_recent > 0 else 0
    
    # Customize the chart
    ax.set_xlabel('Plantae Percentage Range', fontsize=12)
    ax.set_ylabel('Number of Datasets', fontsize=12)
    ax.set_title(f'Distribution of Plantae Content in Recently Modified Datasets\n'
                f'(Last {days} days: {total_with_plants} out of {total_recent} datasets - {plant_percentage:.1f}%)',
                fontsize=14, fontweight='bold', pad=15)
    
    # Add a grid for better readability
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Adjust layout and save
    plt.tight_layout()
    save_figure(fig, 'recent_datasets_plantae_distribution.png')
    
    print(f"\nAnalysis of recently modified datasets:")
    print(f"Total datasets modified in the last {days} days: {total_recent}")
    print(f"Datasets with Plantae records: {total_with_plants} ({plant_percentage:.1f}%)")
    for category, count in range_counts.items():
        print(f"  {category}: {count} datasets")
    
    plt.close()

if __name__ == "__main__":
    print("Filtering datasets at the dataset level (license/type) using random page selection...")
    # Get 500 datasets from each license type
    dataset_keys = get_dataset_keys(COMMERCIAL_LICENSES, page_size=500, datasets_per_license=500)
    print(f"Number of occurrence datasets from random page selection: {len(dataset_keys)}")
    if not dataset_keys:
        print("No datasets found.")
    else:
        # Calculate plant percentages for all datasets
        plantae_percentages = calculate_plantae_percentages(dataset_keys)
        
        print(f"Number of datasets with at least one Plantae record: {len(plantae_percentages)}")
        print("Plantae percentages for all datasets with at least one Plantae record:")
        for k, v in plantae_percentages.items():
            print(f"{k}: {v:.2f}%")
        
        # Create CDF graph of Plantae percentages
        print("\nGenerating cumulative distribution graph of Plantae percentages...")
        create_plantae_percentage_cdf(plantae_percentages)
            
        # Only keep datasets with >50% Plantae
        high_plantae_datasets = [k for k, v in plantae_percentages.items() if v > 50]
        print(f"\nNumber of datasets with >50% Plantae: {len(high_plantae_datasets)}")
        
        if high_plantae_datasets:
            # Take a sample of up to 10 datasets
            sample = random.sample(high_plantae_datasets, min(10, len(high_plantae_datasets)))
            
            # Analyze each dataset individually
            print("\n=== Individual Dataset Analysis ===")
            dataset_stats = []
            for i, selected_dataset in enumerate(sample, 1):
                print(f"\nAnalyzing quality for dataset {i}: {selected_dataset}")
                plant_occurrences = get_plant_occurrences_for_dataset(selected_dataset, n=1000)
                if not plant_occurrences:
                    print("No Plantae records with coordinates found in this dataset.")
                    continue
                num_kingdom_plants = sum(1 for rec in plant_occurrences if rec.get("kingdomKey") == PLANT_TAXON_KEY)
                print(f"Number of records fetched: {len(plant_occurrences)}; Number with kingdomKey==6: {num_kingdom_plants}")
                total, fit, percent = analyze_occurrence_quality(plant_occurrences)
                dataset_stats.append((selected_dataset, total, fit, percent))
            
            # Generate the overall bar graph for quality issue flags
            print("\n=== Cross-Dataset Quality Issues Analysis ===")
            analyze_quality_issues_across_datasets(sample)
            
            # Generate dataset comparison visualizations
            print("\n=== Dataset Comparison Analysis ===")
            dataset_profiles = []
            for dataset_key in sample:
                profile = analyze_dataset_issue_profile(dataset_key)
                if profile:
                    dataset_profiles.append(profile)
            
            if dataset_profiles:
                create_dataset_comparison_chart(dataset_profiles)
            
            # Print summary of dataset quality
            print("\n=== Quality Summary for All Datasets ===")
            for ds_key, total, fit, percent in dataset_stats:
                print(f"Dataset {ds_key}: {percent:.2f}% fit-to-go ({fit}/{total} records)")
        else:
            print("No datasets with >50% Plantae found.")
            
    # Analyze recently modified datasets
    analyze_recently_modified_datasets(days=30, license_types=COMMERCIAL_LICENSES)
    
    # Analyze recently modified datasets
    analyze_recently_modified_datasets(days=30, license_types=COMMERCIAL_LICENSES)