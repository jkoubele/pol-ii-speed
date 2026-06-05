import json
import requests
import time
from urllib.request import urlretrieve
from pathlib import Path


def get_experiment_json(accession: str):
    base_url = "https://www.encodeproject.org"
    url = f"{base_url}/experiments/{accession}/?format=json&frame=embedded"
    response = requests.get(url, headers={"Accept": "application/json"})
    time.sleep(0.1)
    return response.json()


encode_json_folder = Path('/home/jakub/Desktop/pol-ii-speed/epigenetics/ENCODE_jsons')
output_folder = Path('/home/jakub/Desktop/pol-ii-speed/epigenetics/chipseq_data')

for encode_json_file in encode_json_folder.iterdir():
    degron_target = encode_json_file.stem

    with open(encode_json_file) as file:
        encode_json = json.load(file)

    for g in encode_json['@graph']:
        if g['assay_title'] in ('Mint-ChIP-seq', 'TF ChIP-seq'):
            # if g['assay_title']  in ('TF ChIP-seq'):
            target = g['target']
            print(g['assay_title'])
            print(target['label'])
            print(g['accession'])
            experiment_json = get_experiment_json(g['accession'])
            biosamples = [replicate['library']['biosample'] for replicate in experiment_json['replicates']]
            print("Treatments:", biosamples[0].get('treatments'))

            treatment_type = 'degron' if biosamples[0].get('treatments') else 'control'

            # select peaks: prefer IDR over pseudoreplicated
            peak_files = [f for f in experiment_json['files']
                          if f['file_format'] == 'bed'
                          and f['output_type'] in ('IDR thresholded peaks', 'pseudoreplicated peaks')]
            idr_files = [f for f in peak_files if f['output_type'] == 'IDR thresholded peaks']
            peak_file = (idr_files or peak_files or [None])[0]

            bigwig_files = [f for f in experiment_json['files']
                            if f['file_format'] == 'bigWig'
                            and f['output_type'] == 'fold change over control']
            bigwig_file = bigwig_files[0] if bigwig_files else None

            if not peak_file and not bigwig_file:
                continue

            output_subfolder = output_folder / degron_target / target['label'] / treatment_type
            output_subfolder.mkdir(exist_ok=True, parents=True)

            if peak_file:
                dest = output_subfolder / 'peaks.bed.gz'
                if not dest.exists():
                    urlretrieve(peak_file["cloud_metadata"]["url"], dest)

            if bigwig_file:
                dest = output_subfolder / 'fold_change_over_control.bigWig'
                if not dest.exists():
                    urlretrieve(bigwig_file["cloud_metadata"]["url"], dest)
