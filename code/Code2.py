import csv
import re

#Preliminary mechanism search to narrow down manual annotation process

# Define keywords related to antimicrobial peptide mechanisms
KEYWORDS = [
    r'\bpore\b',
    r'\bcarpet\b',
    r'\bmechanism\b',
    r'\bmembrane permeability\b',
    r'\bbarrel-stave\b',
]

# Compile regex patterns
keyword_patterns = [(keyword, re.compile(keyword, re.IGNORECASE)) for keyword in KEYWORDS]

def check_mechanism_keywords(pmc_list_file, input_dir, output_file):
    # Read PMC IDs from the file
    with open(pmc_list_file, 'r') as f:
        pmc_ids = set(line.strip() for line in f)
    
    # Open output file for writing results
    with open(output_file, 'w', newline='', encoding='UTF-8') as f_output:
        writer = csv.writer(f_output, delimiter='\t')
        writer.writerow(['pmcid', 'section', 'paragraph', 'sentence', 'text', 'hit'])  # Header for output
        
        # Counters for statistics
        total_articles = 0
        articles_with_hits = 0
        articles_without_hits = []  # List to track articles with no matches
        
        # Iterate over files in the input directory
        for pmc_id in pmc_ids:
            file_path = f"{input_dir}/{pmc_id}.tsv"  # Assuming each file is named by PMC ID
            article_has_match = False
            total_articles += 1
            
            try:
                with open(file_path, 'r', newline='', encoding='UTF-8') as f_input:
                    reader = csv.reader(f_input, delimiter='\t')
                    
                    # Search each row for keywords
                    for row in reader:
                        section, paragraph, sentence, text, pmcid = row
                        
                        # Check for keyword matches
                        for keyword, pattern in keyword_patterns:
                            if pattern.search(text):
                                # Write matching row with the keyword hit to the output file
                                writer.writerow([pmcid, section, paragraph, sentence, text, keyword])
                                print(f"Match found in article {pmcid}: {text} (Keyword: {keyword})")
                                article_has_match = True  # Mark article as having a match
            
            except FileNotFoundError:
                print(f"File for PMC ID {pmc_id} not found in directory {input_dir}")
            
            # Update count if article has any matches
            if article_has_match:
                articles_with_hits += 1
            else:
                # Add articles without matches to the list
                articles_without_hits.append(pmc_id)
        
        # Calculate percentage of articles with matches
        percentage_with_hits = (articles_with_hits / total_articles) * 100 if total_articles > 0 else 0
        
        # Print the statistics
        print(f"Total articles checked: {total_articles}")
        print(f"Articles with hits: {articles_with_hits} ({percentage_with_hits:.2f}%)")
        print(f"Articles without hits: {len(articles_without_hits)}")
        
        # Print articles with no matches
        print("Articles without matches:")
        for pmc_id in articles_without_hits:
            print(f"PMC ID {pmc_id} has no matching sequences.")

if __name__ == '__main__':
    pmc_list_file = 'PMCarticles_with_matches.txt'  # File containing PMC IDs
    input_dir = 'articles'  # Directory containing article files
    output_file = 'mechanism_keyword_matches.tsv'  # Output file
    
    check_mechanism_keywords(pmc_list_file, input_dir, output_file)
