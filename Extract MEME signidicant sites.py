import json
import pandas as pd

# =====================================
# Load MEME JSON
# =====================================
json_file = r"C:\Users\Dell\Downloads\result.json"

with open(json_file, "r", encoding="utf-8") as f:
    data = json.load(f)

rows = data["MLE"]["content"]["0"]

results = []

for i, row in enumerate(rows, start=1):

    alpha = row[0]
    beta1 = row[1]
    p1 = row[2]
    beta_plus = row[3]
    p_plus = row[4]
    lrt = row[5]
    pvalue = row[6]
    branches = row[7]
    total_branch_length = row[8]
    meme_logl = row[9]
    fel_logl = row[10]
    fel_alpha = row[11]
    fel_beta = row[12]

    if pvalue <= 0.10:

        results.append({
            "Codon": i,
            "Alpha": alpha,
            "Beta+": beta_plus,
            "LRT": lrt,
            "P-value": pvalue,
            "Branches": branches,
            "FEL Alpha": fel_alpha,
            "FEL Beta": fel_beta
        })

df = pd.DataFrame(results)

print(df)
print()
print("Number of significant sites:", len(df))

df.to_excel("MEME_Positive_Selection.xlsx", index=False)
df.to_csv("MEME_Positive_Selection.csv", index=False)

print("\nFiles saved.")