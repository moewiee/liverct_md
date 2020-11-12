import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from tqdm import tqdm
import glob

LOCAL_CSV = "liver_ct_full_export_271020.csv"
GLOBAL_CSV = "liver_ct_global_export_271020.csv"

def parse_xml(label_file):
    """
    Function to parse xml file. Input xml file path, output:
    df: dataframe of local annotation points
    global_labels: global labels in this xml
    studyID: StudyInstanceUID from dicom
    seriesID: SeriesInstanceUID from dicom
    patientpid: PatientPID on label-pacs system
    session_id: Name of annotated doctor
    """
    label_file = label_file
    study_id = label_file.split('/')[-2]
    # Loading label file
    try:
        tree = ET.parse(label_file)
    except FileNotFoundError:
        study_id = label_file.split('/')[-2]
        print(f'Study_id = {study_id} Has no xml file ')
        return None
    
    data = {'studyUid': [],
            'seriesUid': [],
            'imageUid': [],
            'createTimestamp': [],
            'sessionId': [],
            'type': [],
            'annotation': [],
            'name': [],
            'x_pos': [],
            'y_pos': [],
            'z_pos': []}
        
    root = tree.getroot()
    
    global_label = ''
    seriesID = ''
    studyID = ''
    patientpid = ''
    session_id = ''
    
    for p in root.iter('patient'):
        patientpid = p.attrib['pid']
    
    for i, label in enumerate(root.iter('label')):
        if not session_id:
            session_id = label.attrib['sessionId']
        studyID = label.attrib['studyUid']
        if not seriesID:
            seriesID = label.attrib['seriesUid']
        if label.attrib['type'] == 'global':
            for value in label.iter('value'):
                try:
                    if not value.attrib["name"] in global_label:
                        global_label += f'{value.attrib["name"]},'
                except KeyError:
                    continue
                    
        
        for key in data.keys():
            point = label.find('point')
            for j, value in enumerate(point.iter('value')):
                if key == 'name':
                    tag = label.find('tags')
                    label_name = ''
                    for value in tag.iter('value'):
                        label_name += f'{value.attrib["name"]}, '
                    data[key].append(label_name)

                elif key == 'x_pos' or key == 'y_pos' or key == 'z_pos':
                    axis = key.split('_')[0]
                    data[key].append(float(value.attrib[axis]))
                else:
                    data[key].append(label.attrib[key])

    df = pd.DataFrame(data)
    global_labels = global_label.split(',')[:-1]
    accepted_timeStamp = []
    for ts in df["createTimestamp"].unique():
        df_ts = df[df["createTimestamp"] == ts]
        if len(df_ts["z_pos"].unique() == 1):
            accepted_timeStamp.append(ts)
    df = df[df["createTimestamp"].isin(accepted_timeStamp)]
    return df, global_labels, studyID, seriesID, patientpid, session_id

# get list of xml files
# can input many folders at once
xml_folders = [
    "xmls/271020",
]
xml_list = []
for f in xml_folders:
    xml_list += glob.glob(f"{f}/*/*.xml")

# get list of global labels
DISEASE = {
 'No finding': 0,
 'Adenoma': 1,
 'Cystic lesions': 2,
 'FNH': 3,
 'Fatty liver': 4,
 'HCC': 5,
 'Hemangioma': 6,
 'Hepatic atrophy': 7,
 'Hepatic hypertrophy': 8,
 'Intrahepatic Cholangiocarcinoma': 9,
 'Liver abscess': 10,
 'Liver cirrhosis': 11,
 'Metastases': 12,
 'Other': 13,
 'Fasciola hepatica': 14,
 'Perihilar Cholangiocarcinoma': 15,
 'Invalid (Không hợp lệ)': 16,
 'Other liver disease': 17,
 'Biliary dilatation': 18,
 'Multiple hepatic cysts': 19,
 'Choledocholithiasis': 20,
 'Common bile duct stone': 21,
 'Portal hypertension': 22,
 'Simple hepatic cyst': 23,
 'Cholecystitis': 24
}

# get inadherent files
rejected_xml = []
with open("inadherent.txt", 'w') as f:
    # xml without globals
    f.writelines("Error 1: XML without global labels\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if not global_labels:
            rejected_xml.append(xml)
            f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml without local
    f.writelines("Error 2: XML without local labels\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if not len(study_df[study_df["type"] == "local"]):
            if global_labels != ['No finding'] and global_labels != ["Invalid (Không hợp lệ)"] and global_labels != ["Other"]:
                rejected_xml.append(xml)
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml without any labels
    f.writelines("Error 3: XML without any labels\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if not len(study_df[study_df["type"] == "local"]) and not global_labels:
            rejected_xml.append(xml)
            f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml with both NF and other global
    f.writelines("Error 4: XML with both NF and global labels\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if 'No finding' in global_labels:
            if len(global_labels) != 1:
                rejected_xml.append(xml)
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml with global NF but exist local boxes
    f.writelines("Error 5: XML with both NF global label and local bounding boxes\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if global_labels == ['No finding']:
            if len(study_df[study_df["type"] == "local"]):            
                rejected_xml.append(xml)
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml with box only have 1 slice
    f.writelines("Error 6: XML with boxes having only one slice\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        local_df = study_df[study_df["type"] == "local"]
        for ts in local_df['createTimestamp'].unique():
            if len(local_df[local_df['createTimestamp'] == ts]) == 4:
                rejected_xml.append(xml)
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml without seriesid
    f.writelines("Error 7: XML without seriesuid\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if not seriesID:
            if global_labels != ["No finding"] and global_labels != ["Invalid (Không hợp lệ)"] and global_labels != ["Other"]:
                rejected_xml.append(xml)
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')
    # xml with 2d labels but no 3d labels
    f.writelines("Error 8: XML without any labels but 2D bounding boxes\n")
    for xml in tqdm(xml_list):
        study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
        if not len(study_df[study_df["annotation"].isin(['20'])]):
            if len(study_df[study_df["annotation"].isin(['Rectangle'])]):
                f.writelines(session_id + ' ' * (10-len(session_id)) + '| ' + xml.split('/')[2] + ' ' * (60 - len(xml.split('/')[2])) + '| ' + patientpid + '\n')

# remove xmls that are not adhere to labelling rules before generate csv file
xml_list = [x for x in xml_list if x not in list(set(rejected_xml))]

# generate local csv
df = pd.DataFrame()
for xml in tqdm(xml_list):
    sub_data = []
    study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
    df = pd.concat([df, study_df])
df = df[df["type"] == "local"]
df = df[df["annotation"] == "20"]
df.to_csv(LOCAL_CSV, index=False)

# generate global csv
data = []
for xml in tqdm(xml_list):
    sub_data = []
    study_df, global_labels, studyID, seriesID, patientpid, session_id = parse_xml(xml)
    if not studyID:
        sub_data.append(' ')
    else:
        sub_data.append(studyID)
    if not seriesID:
        sub_data.append(' ')
    else:
        sub_data.append(seriesID)
    for _ in range(len(DISEASE)):
        sub_data.append(0)
    for lb in global_labels:
        sub_data[DISEASE[lb]+2] = 1
    if not sub_data[DISEASE["Invalid (Không hợp lệ)"]+2]:
        if not (sub_data[DISEASE["Other"]+2] and np.sum(sub_data[2:]) == 1 and not seriesID):
            data.append(sub_data)    
df = pd.DataFrame(data=data, columns=['StudyID', 'SeriesID'] + list(DISEASE.keys()))
df.to_csv(GLOBAL_CSV, index=False)