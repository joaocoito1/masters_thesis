
"""
Builds an Excel for manually classifying HGT candidates. Contains a
Groups sheet (one row per tree, with a Clear HGT/Clear Not HGT/Unclear verdict) 
and a Candidates sheet (one row per protein, verdict auto-filled from its group, 
with manual override).
"""


import os
import pandas as pd
import datetime

import openpyxl
from openpyxl.styles import PatternFill
from openpyxl.utils import get_column_letter
from openpyxl.formatting.rule import CellIsRule
from openpyxl.worksheet.table import Table, TableStyleInfo
from openpyxl.worksheet.datavalidation import DataValidation


def generate_review_excel(strains_df, log, AlienIndex_path, hgt_groups_path, groups_wo_trees_path, hgt_review_path):

    # Extract strain names to be analyzed
    strain_names = strains_df["strain_name"].tolist()


    ##### Load all Alien Index results into one DataFrame #####
    AI_dfs = []
    for strain in strain_names:
        # Build input file path and check that it exists
        AI_file = os.path.join(AlienIndex_path, f"{strain}_AlienIndex_results.csv")

        AI_result = pd.read_csv(AI_file)
        AI_result["strain"] = strain
        AI_result = AI_result[["strain"] + [c for c in AI_result.columns if c != "strain"]] # Put "strain" in 1st column
        AI_dfs.append(AI_result)

    # Join all results
    AlienIndex_results = pd.concat(AI_dfs, ignore_index=True)


    ##### Add Group information #####

    # Get groups
    with open(hgt_groups_path) as f:
        group_to_query = {}
        for line in f.readlines():
            group_number, group = line.strip().split(": ")
            group = set(group.split())
            group_to_query[group_number] = group

    # Add group column
    query_to_group = {q: g_number for g_number, members in group_to_query.items() for q in members}
    AlienIndex_results["Group"] = AlienIndex_results["query"].map(query_to_group)

    # Keep only the ones that belong to a group
    candidates_table = AlienIndex_results.dropna(subset=["Group"]).copy()

    # Count candidates per group
    group_table = (candidates_table
                    .groupby("Group")["query"]
                    .nunique()
                    .reset_index(name="Group size"))


    ##### Build excel #####

    # Group sheet
    group_table["Group verdict"] = ""
    group_table["Donor"] = ""
    group_table["Notes"] = ""
    group_table = group_table.sort_values("Group size", ascending=False).reset_index(drop=True)

    # Mark groups without a tree as Unclear, with a note
    with open(groups_wo_trees_path) as f:
        missing_groups = {line.strip() for line in f if line.strip()}
    no_tree_group = group_table["Group"].isin(missing_groups)
    group_table.loc[no_tree_group, "Group verdict"] = "Unclear"
    group_table.loc[no_tree_group, "Notes"] = "No Tree"

    # Candidate sheet
    candidates_table["Group verdict"] = ""       # will be a formula
    candidates_table["Individual verdict"] = ""  # will be dropdown
    candidates_table["Final verdict"] = ""       # will be a formula
    
    candidates_table["Group donor"] = ""         # will be a formula
    candidates_table["Individual donor"] = ""    # will be dropdown
    candidates_table["Final donor"] = ""         # will be a formula

    candidates_table["Notes"] = ""               # Free text


    # Writte to excel
    with pd.ExcelWriter(hgt_review_path, engine="openpyxl") as writer:
        group_table.to_excel(writer, sheet_name="Groups", index=False)
        candidates_table.to_excel(writer, sheet_name="Candidates", index=False)


    ##### Add formulas #####

    # Load excel
    wb = openpyxl.load_workbook(hgt_review_path)
    sheet_g = wb["Groups"]
    sheet_c = wb["Candidates"]

    # Column positions
    col_g = {name: group_table.columns.get_loc(name) + 1 for name in group_table.columns}
    col_c = {name: candidates_table.columns.get_loc(name) + 1 for name in candidates_table.columns}
    last_row_g = sheet_g.max_row
    last_row_c = sheet_c.max_row

    group_letter = get_column_letter(col_c["Group"])
    gv_letter = get_column_letter(col_c["Group verdict"])
    iv_letter = get_column_letter(col_c["Individual verdict"])
    gd_letter = get_column_letter(col_c["Group donor"])
    id_letter = get_column_letter(col_c["Individual donor"])

    groups_lookup_range_verdict = f"Groups!$A$2:${get_column_letter(col_g['Group verdict'])}${last_row_g}"
    groups_lookup_range_donor = f"Groups!$A$2:${get_column_letter(col_g['Donor'])}${last_row_g}"

    for r in range(2, last_row_c + 1):
        sheet_c.cell(row=r, column=col_c["Group verdict"]).value = (
            f'=IFERROR(VLOOKUP({group_letter}{r},{groups_lookup_range_verdict},{col_g["Group verdict"]},FALSE),"")')

        sheet_c.cell(row=r, column=col_c["Final verdict"]).value = (
            f'=IF({iv_letter}{r}<>"",{iv_letter}{r},'
            f'IF({gv_letter}{r}="Clear HGT","HGT",'
            f'IF({gv_letter}{r}="Clear Not HGT","Not HGT","PENDING REVIEW")))')

        sheet_c.cell(row=r, column=col_c["Group donor"]).value = (
            f'=IFERROR(VLOOKUP({group_letter}{r},{groups_lookup_range_donor},{col_g["Donor"]},FALSE),"")')

        sheet_c.cell(row=r, column=col_c["Final donor"]).value = (
            f'=IF({id_letter}{r}<>"",{id_letter}{r},'
            f'IF({gd_letter}{r}<>"",{gd_letter}{r},"PENDING REVIEW"))')

    ##### Add conditional formatting #####

    GREEN, RED, YELLOW = "C6EFCE", "FFC7CE", "FFEB9C"
    BLUE, ORANGE, PURPLE = "B4C7E7", "F8CBAD", "D9D2E9"

    # Groups sheet: color Group verdict
    gv_g_letter = get_column_letter(col_g["Group verdict"])
    gv_g_range = f"{gv_g_letter}2:{gv_g_letter}{last_row_g}"
    sheet_g.conditional_formatting.add(gv_g_range, CellIsRule(operator="equal", formula=['"Clear HGT"'], fill=PatternFill(bgColor=GREEN)))
    sheet_g.conditional_formatting.add(gv_g_range, CellIsRule(operator="equal", formula=['"Clear Not HGT"'], fill=PatternFill(bgColor=RED)))
    sheet_g.conditional_formatting.add(gv_g_range, CellIsRule(operator="equal", formula=['"Unclear"'], fill=PatternFill(bgColor=YELLOW)))

    # Groups sheet: color Donor
    donor_g_letter = get_column_letter(col_g["Donor"])
    donor_g_range = f"{donor_g_letter}2:{donor_g_letter}{last_row_g}"
    sheet_g.conditional_formatting.add(donor_g_range, CellIsRule(operator="equal", formula=['"Fungi"'], fill=PatternFill(bgColor=BLUE)))
    sheet_g.conditional_formatting.add(donor_g_range, CellIsRule(operator="equal", formula=['"Bacteria"'], fill=PatternFill(bgColor=ORANGE)))
    sheet_g.conditional_formatting.add(donor_g_range, CellIsRule(operator="equal", formula=['"Mixed"'], fill=PatternFill(bgColor=PURPLE)))

    # Candidates sheet: color Group verdict (echoed from Groups sheet)
    gv_c_range = f"{gv_letter}2:{gv_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(gv_c_range, CellIsRule(operator="equal", formula=['"Clear HGT"'], fill=PatternFill(bgColor=GREEN)))
    sheet_c.conditional_formatting.add(gv_c_range, CellIsRule(operator="equal", formula=['"Clear Not HGT"'], fill=PatternFill(bgColor=RED)))
    sheet_c.conditional_formatting.add(gv_c_range, CellIsRule(operator="equal", formula=['"Unclear"'], fill=PatternFill(bgColor=YELLOW)))

    # Candidates sheet: color Individual verdict
    iv_range = f"{iv_letter}2:{iv_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(iv_range, CellIsRule(operator="equal", formula=['"HGT"'], fill=PatternFill(bgColor=GREEN)))
    sheet_c.conditional_formatting.add(iv_range, CellIsRule(operator="equal", formula=['"Not HGT"'], fill=PatternFill(bgColor=RED)))

    # Candidates sheet: color Final verdict
    fv_letter = get_column_letter(col_c["Final verdict"])
    fv_range = f"{fv_letter}2:{fv_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(fv_range, CellIsRule(operator="equal", formula=['"HGT"'], fill=PatternFill(bgColor=GREEN)))
    sheet_c.conditional_formatting.add(fv_range, CellIsRule(operator="equal", formula=['"Not HGT"'], fill=PatternFill(bgColor=RED)))
    sheet_c.conditional_formatting.add(fv_range, CellIsRule(operator="equal", formula=['"PENDING REVIEW"'], fill=PatternFill(bgColor=YELLOW)))

    # Candidates sheet: color Group donor (echoed from Groups sheet)
    gd_range = f"{gd_letter}2:{gd_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(gd_range, CellIsRule(operator="equal", formula=['"Fungi"'], fill=PatternFill(bgColor=BLUE)))
    sheet_c.conditional_formatting.add(gd_range, CellIsRule(operator="equal", formula=['"Bacteria"'], fill=PatternFill(bgColor=ORANGE)))
    sheet_c.conditional_formatting.add(gd_range, CellIsRule(operator="equal", formula=['"Mixed"'], fill=PatternFill(bgColor=PURPLE)))

    # Candidates sheet: color Individual donor
    id_range = f"{id_letter}2:{id_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(id_range, CellIsRule(operator="equal", formula=['"Fungi"'], fill=PatternFill(bgColor=BLUE)))
    sheet_c.conditional_formatting.add(id_range, CellIsRule(operator="equal", formula=['"Bacteria"'], fill=PatternFill(bgColor=ORANGE)))
    sheet_c.conditional_formatting.add(id_range, CellIsRule(operator="equal", formula=['"Mixed"'], fill=PatternFill(bgColor=PURPLE)))

    # Candidates sheet: color Final donor
    fd_letter = get_column_letter(col_c["Final donor"])
    fd_range = f"{fd_letter}2:{fd_letter}{last_row_c}"
    sheet_c.conditional_formatting.add(fd_range, CellIsRule(operator="equal", formula=['"Fungi"'], fill=PatternFill(bgColor=BLUE)))
    sheet_c.conditional_formatting.add(fd_range, CellIsRule(operator="equal", formula=['"Bacteria"'], fill=PatternFill(bgColor=ORANGE)))
    sheet_c.conditional_formatting.add(fd_range, CellIsRule(operator="equal", formula=['"Mixed"'], fill=PatternFill(bgColor=PURPLE)))
    sheet_c.conditional_formatting.add(fd_range, CellIsRule(operator="equal", formula=['"PENDING REVIEW"'], fill=PatternFill(bgColor=YELLOW)))

    ##### Add table #####

    def add_table(ws, name, n_rows, n_cols):
        last_col_letter = get_column_letter(n_cols)
        table_range = f"A1:{last_col_letter}{n_rows}"
        table = Table(displayName=name, ref=table_range)
        table.tableStyleInfo = TableStyleInfo(
            name="TableStyleMedium2", showRowStripes=True
        )
        ws.add_table(table)

    add_table(sheet_g, "GroupsTable", last_row_g, len(group_table.columns))
    add_table(sheet_c, "CandidatesTable", last_row_c, len(candidates_table.columns))

    ##### Add dropdown menus #####

    # Groups sheet
    dv_group = DataValidation(type="list", formula1='"Clear HGT,Clear Not HGT,Unclear"', allow_blank=True)
    sheet_g.add_data_validation(dv_group)
    gv_g_letter = get_column_letter(col_g["Group verdict"])
    dv_group.add(f"{gv_g_letter}2:{gv_g_letter}{last_row_g}")

    dv_donor = DataValidation(type="list", formula1='"Fungi,Bacteria,Mixed"', allow_blank=True)
    sheet_g.add_data_validation(dv_donor)
    donor_g_letter = get_column_letter(col_g["Donor"])
    dv_donor.add(f"{donor_g_letter}2:{donor_g_letter}{last_row_g}")

    # Candidates sheet
    dv_individual_verdict = DataValidation(type="list", formula1='"HGT,Not HGT"', allow_blank=True)
    sheet_c.add_data_validation(dv_individual_verdict)
    dv_individual_verdict.add(f"{iv_letter}2:{iv_letter}{last_row_c}")

    dv_individual_donor = DataValidation(type="list", formula1='"Fungi,Bacteria,Mixed"', allow_blank=True)
    sheet_c.add_data_validation(dv_individual_donor)
    dv_individual_donor.add(f"{id_letter}2:{id_letter}{last_row_c}")

    ##### Extra #####

    # Hide specific columns
    for col_name in ["strain","outgroup_norm_bitscore","ingroup_norm_bitscore","self_bitscore","outgroup_subject","ingroup_subject"]:
        col_letter = get_column_letter(col_c[col_name])
        sheet_c.column_dimensions[col_letter].hidden = True

    # Define decimal places shown for AlienIndex
    col_letter = get_column_letter(col_c["AI"])
    for cell in sheet_c[col_letter][1:]:    # [1:] -> skip header
        cell.number_format = "0.00"

    wb.save(hgt_review_path)
    log(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} : {os.path.basename(hgt_review_path)} created.")
    
