import csv
import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt

#Extract peptide sequences from articles

# Define regex patterns for amino acid sequences
REGEXES = [
    r'[AC-IK-NP-TVWY]{10,60}'  # Matches amino acid sequences between 10 and 60 residues
]

# Mapping specific phrases to general categories
SECTION_CATEGORIES = {
    "materials and methods": "Materials and Methods",
    "methods": "Materials and Methods",
    "results": "Results",
    "discussion": "Discussion",
    "introduction": "Introduction",
    "abstract": "Abstract"
}

def categorize_section(section):
    section_lower = section.lower()
    for key, category in SECTION_CATEGORIES.items():
        if key in section_lower:
            return category
    return "Other"

def is_dna_sequence(sequence):
    """Checks if a sequence contains only DNA bases (A, G, C, T)."""
    return all(char in "AGCT" for char in sequence)

def regex_search():
    regexes_compiled = [re.compile(regex) for regex in REGEXES]
    dir_export = 'articles'
    
    # Output files
    name_output_matches = 'amino_acid_matches.tsv'
    name_output_summary = 'summary_statistics.txt'
    
    # Open output file for matches
    with open(name_output_matches, 'w', newline='', encoding='UTF-8') as f_output_matches:
        filewriter = csv.writer(f_output_matches, delimiter='\t')
        header = ['section', 'paragraph', 'sentence', 'text', 'pmcid', 'matched_sequence']
        filewriter.writerow(header)
        
        # Counters and data structures for summary statistics
        total_articles = 0
        articles_with_matches = set()  # To store unique articles with matches
        section_counts = defaultdict(int)  # Counts of sections with matches
        pmcid_list = []  # List of PMC IDs with relevant matches
        
        # Iterate over files in the directory
        for obj in os.listdir(dir_export):
            path_input = os.path.join(dir_export, obj)
            with open(path_input, 'r', newline='', encoding='UTF-8') as f_input:
                total_articles += 1
                filereader = csv.reader(f_input, delimiter='\t')
                has_match = False
                
                # Process each row in the file
                for row in filereader:
                    section, paragraph, sentence, text, pmcid = row
                    
                    # Generalize the section category
                    general_section = categorize_section(section)
                    
                    # Check for regex matches
                    for re_comp in regexes_compiled:
                        matches = re.findall(re_comp, text)
                        for sequence in matches:
                            # Filter out DNA sequences (A, G, C, T only)
                            if not is_dna_sequence(sequence):
                                # Write each matched sequence as a new row
                                filewriter.writerow([section, paragraph, sentence, text, pmcid, sequence])  
                                
                                # Update summary data
                                articles_with_matches.add(pmcid)
                                pmcid_list.append(pmcid)
                                section_counts[general_section] += 1
                                has_match = True
                
                if has_match:
                    articles_with_matches.add(pmcid)
        
        # Summary statistics
        matched_articles_count = len(articles_with_matches)
        percentage_with_matches = (matched_articles_count / total_articles) * 100
        
        # Write summary statistics to a separate file
        with open(name_output_summary, 'w', encoding='UTF-8') as f_output_summary:
            f_output_summary.write(f"Total articles processed: {total_articles}\n")
            f_output_summary.write(f"Articles with matches: {matched_articles_count} ({percentage_with_matches:.2f}%)\n")
            f_output_summary.write("PMC IDs with matches:\n")
            for pmcid in articles_with_matches:
                f_output_summary.write(f"{pmcid}\n")
            f_output_summary.write("\nFrequency of matches by section:\n")
            for section, count in section_counts.items():
                f_output_summary.write(f"{section}: {count}\n")
    
    # Plotting histogram of section frequencies
    sections = list(section_counts.keys())
    frequencies = list(section_counts.values())
    
    plt.figure(figsize=(10, 6))
    plt.bar(sections, frequencies)
    plt.xlabel("Section")
    plt.ylabel("Frequency of Matches")
    plt.title("Frequency of Matches by Section in Articles")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("section_frequency_histogram.png")
    plt.show()

if __name__ == '__main__':
    regex_search()
