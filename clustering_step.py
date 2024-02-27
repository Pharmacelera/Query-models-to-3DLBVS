import numpy as np
from rdkit import Chem
import argparse
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

def dihedrals_angles(mol,circular = True):

    raw_rot_bonds =  mol.GetSubstructMatches(Chem.MolFromSmarts("[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]"))
    bonds = []
    rot_bonds = []

    for k,i,j,l in raw_rot_bonds:
        if (i,j) not in bonds:
            bonds.append((i,j))
            rot_bonds.append((k,i,j,l))
    thetas = []
    dih_points = []
    if len(rot_bonds) == 0:
        print('Not rotatable bonds found. Exiting clustering protocol ...')
        sys.exit()

    for k,i,j,l in rot_bonds:
        if (k+1,i+1,j+1,l+1) not in dih_points:
            dih_points.append((k+1,i+1,j+1,l+1))

        dih = Chem.rdMolTransforms.GetDihedralRad(mol.GetConformers()[0], k,i,j,l )

        if circular:
            thetas.append(np.cos(dih))
            thetas.append(np.sin(dih))
        else:
            thetas.append(dih)

    return thetas,dih_points

def get_dPCA_mat(mols):

    dihedrals = []
    dih_mat = []
    for mol in mols:
        if mol != None:
            dihedrals,_ = dihedrals_angles(mol,circular=True)
            dih_mat.append(dihedrals)
    dih_mat = np.array(dih_mat)
    num_pc = len(dihedrals)
    X = StandardScaler().fit_transform(dih_mat) #normalizing the features
    pca = PCA(n_components=0.75)
    X_pca = pca.fit_transform(X)
    print('N confs: {} and N features: {}'.format(X_pca.shape[0],X_pca.shape[1]))
    return X_pca

def ClusterIndices(clustNum, labels_array):
    return np.where(labels_array == clustNum)[0]

def select_lowest_energy_centroid(model,energy,nclusters):
    min_clusters = []
    for i in range(nclusters):
        idx_cluster = ClusterIndices(i,model.labels_)
        cluster_energy = [(energy[i],i) for i in idx_cluster]
        en_min = 10e6
        for en,id_conf in cluster_energy:
            if en < en_min:
                en_min = en
                min_cluster_idx = id_conf
        min_clusters.append(min_cluster_idx)
    return min_clusters

def get_representatives(dPCA_mat,energy):
    model = KMeans(n_clusters=4)
    model.fit(dPCA_mat)
    # Apply vector quantization to find closest point or "conformer" to the centroid of a particular cluster
    closest = select_lowest_energy_centroid(model,energy,4)
    return closest

def best_energy_id(mols):
    en_min = 10e6
    for mol in mols:
        if mol.HasProp('PHA_CONF_ID'):
            id_mol = mol.GetProp('PHA_CONF_ID')
            # energy provided by PharmScreen
            en = float(mol.GetProp('PHA_CONF_ENERGY').split(' ')[0])
            if en < en_min:
                en_min = en
                id_mol_min = id_mol
                low_en_mol = mol
        else:
            print(" The input sdf has not include 'PHA_CONF_ID', be sure that you are using the proper Pharmscreen output file")
            quit()

    return id_mol_min,low_en_mol

def cluster_conformers(input,output):
    mols = Chem.SDMolSupplier(input,removeHs=False)
    ensemble_PHA_ids = []
    w = Chem.SDWriter(output)
    best_conf_id,le_mol = best_energy_id(mols)
    w.write(le_mol)
    print('Lowest energy conformer PHA id: ',best_conf_id)
    print("\n")
    mols_filtered = []
    for mol in enumerate(mols):
        mols_filtered.append(mol)

    mols2 = [mol for mol in mols_filtered if mol.GetProp('PHA_CONF_ID')!= best_conf_id]

    energy = [float(mol.GetProp('PHA_CONF_ENERGY').split(' ')[0]) for mol in mols2]

    dPCA_mat = get_dPCA_mat(mols2)
    representatives = get_representatives(dPCA_mat,energy)
    for id in representatives:
        ensemble_PHA_ids.append(mols2[id].GetProp('PHA_CONF_ID'))
        w.write(mols2[id])
    return ensemble_PHA_ids

def main():
    parser = argparse.ArgumentParser(description='Script to cluster conformers from multi-coformer pharmscreen sdf')
    parser.add_argument('-i', required=True, help='sdf input file')
    parser.add_argument('-o', required=True, help='sdf output file')

    args = parser.parse_args()
    medoids = cluster_conformers(args.i,args.o)
    print("Representatives ids: ",medoids)

if __name__ == '__main__':
    main()
