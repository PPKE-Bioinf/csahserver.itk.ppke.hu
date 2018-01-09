#!/usr/bin/python3

import ftputil
import os
from datetime import datetime

import yaml

dir_name = os.path.dirname(os.path.realpath(__file__))

with open(dir_name + "/config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

lastversion_path = dir_name + "/lastversion.txt"

with ftputil.FTPHost("ftp.ebi.ac.uk", "anonymous", "anonymous") as ftp_host:
    ftp_path_trembl = "pub/databases/uniprot/knowledgebase/uniprot_trembl.fasta.gz"
    ftp_path_sprot = "pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
    current_timestamp = ftp_host.path.getmtime(ftp_path_trembl)

    try:
        with open(lastversion_path) as fp:  
            saved_timestamp = float(fp.readline())
            print(saved_timestamp)
            print(current_timestamp)
            if saved_timestamp == current_timestamp:
                print("No new upload found.")
                raise SystemExit(0)
    except FileNotFoundError:
        print("Initial startup.")

    print("New upload found.")

    with open(lastversion_path, "w") as lastversion_file:
        lastversion_file.write(str(current_timestamp))

    ftp_host.download(ftp_path_trembl, str(current_timestamp) + "_uniprot_trembl.fasta.gz")
    ftp_host.download(ftp_path_sprot, str(current_timestamp) + "_uniprot_sprot.fasta.gz")


sent_from = "bioinf.notification@gmail.com@gmail.com"
to = cfg['addresses'] 
subject = "OMG Super Important Message"
body = "Subject: CSAH server notification\n\nNew database is available!"

gmail_user = "bioinf.notification@gmail.com"
gmail_password = cfg['password']


import smtplib

try:  
    server = smtplib.SMTP("smtp.gmail.com", 587)
    server.ehlo()
    server.starttls()
    server.ehlo()
    server.login(gmail_user, gmail_password)
    server.sendmail(sent_from, to, body)
    server.close()
    print("Notification e-mail sent.")

except Exception as e:  
    print("[ERROR] e-mail was not sent", e)