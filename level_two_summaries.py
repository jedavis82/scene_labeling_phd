"""
Compute the general domain and person domain level two summaries using the set of level one summaries as input.
Each of the level two summary types will be stored in individual CSV files
"""
import pandas as pd
from level_two_utils.general_rules import GeneralRules
from level_two_utils.person_rules import PersonRules
from tqdm import tqdm
import json


def get_consensus_angle(f0, f2, hyb):
    if f0 == f2:
        return f0
    if f0 == hyb:
        return hyb
    if f2 == hyb:
        return hyb
    else:
        return hyb


def convert_meta_labels(labels):
    for l in labels:
        if l == 'soccer_ball' or l == 'tennis_ball' or l == 'baseball':
            return l
    else:
        return labels[0]


def compute_general_summaries(df, output_file):
    general_rules = GeneralRules()
    results = []
    for idx, row in tqdm(df.iterrows(), total=len(df)):
        key = row['key']
        rel_path = row['relative_path']
        img_name = row['img_name']
        arg_label = row['arg_label']
        arg_bounding_box = row['arg_bounding_box']
        ref_label = row['ref_label']
        ref_bounding_box = row['ref_bounding_box']
        proximity = row['proximity']
        overlap = row['overlap']
        f0 = row['f0']
        f2 = row['f2']
        hybrid = row['hybrid']
        l1_summary = row['level_one_summary']
        interacting_label = general_rules.compute_interactions(proximity, overlap)
        l2_summary = f'{arg_label} {interacting_label} {ref_label}'
        results.append({
            'key': key,
            'relative_path': rel_path,
            'img_name': img_name,
            'arg_label': arg_label,
            'arg_bounding_box': arg_bounding_box,
            'ref_label': ref_label,
            'ref_bounding_box': ref_bounding_box,
            'proximity': proximity,
            'overlap': overlap,
            'f0': f0,
            'f2': f2,
            'hybrid': hybrid,
            'level_one_summary': l1_summary,
            'general_level_two_summary': l2_summary
        })
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, encoding='utf-8', header=True, index=False)
    return results_df


def compute_person_summaries(df, ontology_df, meta_df, output_file):
    person_rules = PersonRules(ontology_df=ontology_df, show_sim_result=False)
    results = []
    for idx, row in tqdm(df.iterrows(), total=len(df)):
        arg_label = row['arg_label']
        if 'person' not in arg_label:
            continue  # No person detection, cannot compute a person domain level two summary
        key = row['key']
        rel_path = row['relative_path']
        img_name = row['img_name']
        arg_bounding_box = row['arg_bounding_box']
        ref_label = row['ref_label']
        ref_bounding_box = row['ref_bounding_box']
        proximity = row['proximity']
        overlap = row['overlap']
        f0 = row['f0']
        f2 = row['f2']
        hybrid = row['hybrid']
        l1_summary = row['level_one_summary']

        # Compute the level two summary for this result
        sr_angle = get_consensus_angle(f0, f2, hybrid)
        meta = json.loads(meta_df.loc[meta_df['relative_path'] == rel_path]['labels'].iloc[0])
        meta_label = convert_meta_labels(meta)
        l2_interaction = person_rules.compute_interactions(ref_label, proximity, overlap, sr_angle, meta_label)
        if l2_interaction is not None:
            l2_summary = f'{arg_label} {l2_interaction} {ref_label}'
        else:
            l2_summary = 'Negative Interaction'

        results.append({
            'key': key,
            'relative_path': rel_path,
            'img_name': img_name,
            'arg_label': arg_label,
            'arg_bounding_box': arg_bounding_box,
            'ref_label': ref_label,
            'ref_bounding_box': ref_bounding_box,
            'proximity': proximity,
            'overlap': overlap,
            'f0': f0,
            'f2': f2,
            'hybrid': hybrid,
            'level_one_summary': l1_summary,
            'person_level_two_summary': l2_summary
        })
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, encoding='utf-8', header=True, index=False)
    return results_df
