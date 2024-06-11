import requests
import os, sys

count = 0
size = len(Protein_id)
print(size)

for id in Protein_id:
    url = "https://www.uniprot.org/uniprot/" + id + ".fasta"

    try:
        r = requests.get(url,timeout=10)
        r.raise_for_status()
    except requests.exceptions.HTTPError as errh:
        print ("Http Error:",errh)
    except requests.exceptions.ConnectionError as errc:
        print ("Error Connecting:",errc)
    except requests.exceptions.Timeout as errt:
        print ("Timeout Error:",errt)
    except requests.exceptions.RequestException as err:
        print ("OOps: Something Else",err)


    dir = os.getcwd() + "/donwload.fasta"

    if r.status_code == 200:
        count += 1
        print(id + " downloaded Successfully")
        with open(dir,'ab+') as f:
            f.write(r.content)
