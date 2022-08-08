import pandas as pd
import numpy as np
from numpy import dot
from numpy.linalg import norm
from scipy.spatial.distance import jaccard
from tqdm import tqdm, tqdm_notebook
from scipy import stats
from scipy.stats import linregress
from scipy.stats import pearsonr

def get_jaccard_correlation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Jaccard similarity for each drug."""
    correlation_matrix = []
    
    for drug in target_matrix.index.values:
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.loc[drug, :]
        dissimilarity_score = jaccard(target_vector, transcriptomic_vector)
        similarity = 1 - dissimilarity_score
        correlation_matrix.append(
            {
                'drug': drug,
                'jaccard_correlation': similarity
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_cosine_correlation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Cosine similarity for each drug."""
    correlation_matrix = []
    
    for drug in target_matrix.index.values:
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.loc[drug, :]
        corr_val = dot(target_vector, transcriptomic_vector)/(norm(target_vector)*norm(transcriptomic_vector))
        correlation_matrix.append(
            {
                'drug': drug,
                'cosine_correlation': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_pearson_corr(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Pearson correlation for each drug."""
    correlation_matrix = []
    
    for drug in target_matrix.index.values:
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.loc[drug, :]
        corr_val, p_val = pearsonr(target_vector, transcriptomic_vector)
        correlation_matrix.append(
            {
                'drug': drug,
                'pearson_correlation': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_pearson_corr_permutation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Pearson correlation for each drug."""
    correlation_matrix = []
    
    for i,drug in enumerate(target_matrix.index.values):
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.iloc[i, :]
        #print(transcriptomic_vector)
        corr_val, p_val = pearsonr(target_vector, transcriptomic_vector)
        correlation_matrix.append(
            {
                'drug': drug,
                'pearson_correlation': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_cosine_corr_permutation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Cosine correlation for each drug."""
    correlation_matrix = []
    
    for i,drug in enumerate(target_matrix.index.values):
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.iloc[i, :]
        corr_val = dot(target_vector, transcriptomic_vector)/(norm(target_vector)*norm(transcriptomic_vector))
        correlation_matrix.append(
            {
                'drug': drug,
                'cosine_correlation': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_jaccard_corr_permutation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Jaccard correlation for each drug."""
    correlation_matrix = []
    
    for i,drug in enumerate(target_matrix.index.values):
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.iloc[i, :]
        #print(transcriptomic_vector)
        dissimilarity_score = jaccard(target_vector, transcriptomic_vector)
        similarity = 1 - dissimilarity_score
        correlation_matrix.append(
            {
                'drug': drug,
                'jaccard_correlation': similarity
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_similarity(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get Pearson correlation for each drug."""
    correlation_matrix = []
    
    for drug in tqdm(target_matrix.index.values):
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.loc[drug, :]
        corr_val = abs(sum(abs(np.array(target_vector) - np.array(transcriptomic_vector))\
                       - np.ones(len(target_vector)))/len(target_vector))
        correlation_matrix.append(
            {
                'drug': drug,
                'similarity': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_similarity_permutation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get similarity for each drug."""
    correlation_matrix = []
    
    for i,drug in enumerate(target_matrix.index.values):
        target_vector = target_matrix.loc[drug, :]
        transcriptomic_vector = transcriptomic_matrix.iloc[i, :]
        #print(transcriptomic_vector)
        corr_val= abs(sum(abs(np.array(target_vector) - np.array(transcriptomic_vector))\
                      - np.ones(len(target_vector)))/len(target_vector))
        correlation_matrix.append(
            {
                'drug': drug,
                'similarity': corr_val
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

from itertools import groupby

def get_nonzero_similarity(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get nonzero_similarity for each drug."""
    correlation_matrix = []
    
    for drug in tqdm(target_matrix.index.values):
        old_target_vector = np.array(target_matrix.loc[drug, :])
        old_transcriptomic_vector = np.array(transcriptomic_matrix.loc[drug, :])
        nonzero_index_1 = np.nonzero(old_target_vector)[0]
        shared_nonzero_index = list(nonzero_index_1)

        #Get values in list that correspond to nonzero indexes
        nonzero_list_1 = [old_target_vector[i] for i in shared_nonzero_index]
        nonzero_list_2 = [old_transcriptomic_vector[i] for i in shared_nonzero_index]

        target_vector = np.array(nonzero_list_1)
        transcriptomic_vector = np.array(nonzero_list_2)
        if list(target_vector):
            corr_val = abs(sum(abs(np.sign(np.array(target_vector) - np.array(transcriptomic_vector)))\
                           - np.ones(len(shared_nonzero_index)))/len(shared_nonzero_index))
        else:
            corr_val = -1
        cosine_corr = dot(target_vector, transcriptomic_vector)/(norm(target_vector)*norm(transcriptomic_vector))
        jaccard_corr = 1 - jaccard(target_vector, transcriptomic_vector)
     
        correlation_matrix.append(
            {
                'drug': drug,
                'nonzero_similarity': corr_val,
                'cosine' : cosine_corr, 
                'jaccard' : jaccard_corr
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df

def get_nonzero_similarity_permutation(
    target_matrix: pd.DataFrame,
    transcriptomic_matrix: pd.DataFrame,
):
    """Get nonzero similarity for each drug."""
    correlation_matrix = []
    
    for i,drug in enumerate(target_matrix.index.values):
        old_target_vector = np.array(target_matrix.loc[drug, :])
        old_transcriptomic_vector = np.array(transcriptomic_matrix.iloc[i, :])
        nonzero_index_1 = np.nonzero(old_target_vector)[0]
        nonzero_index_2 = np.nonzero(old_transcriptomic_vector)[0]
        shared_nonzero_index = list(nonzero_index_1)
        groups = groupby(enumerate(shared_nonzero_index), lambda x: x[1]-x[0])
        groups = [list(y) for x,y in groups]
        
        #Get values in list that correspond to nonzero indexes
        nonzero_list_1 = [old_target_vector[i] for i in shared_nonzero_index]
        nonzero_list_2 = [old_transcriptomic_vector[i] for i in shared_nonzero_index]

        target_vector = np.array(nonzero_list_1)
        transcriptomic_vector = np.array(nonzero_list_2)
        if list(target_vector):
            corr_val = abs(sum(abs(np.sign(np.array(target_vector) - np.array(transcriptomic_vector)))\
                           - np.ones(len(shared_nonzero_index)))/len(shared_nonzero_index))
        else:
            corr_val = -1
        cosine_corr = dot(target_vector, transcriptomic_vector)/(norm(target_vector)*norm(transcriptomic_vector))
        jaccard_corr = 1 - jaccard(target_vector, transcriptomic_vector)
        correlation_matrix.append(
            {
                'drug': drug,
                'nonzero_similarity': corr_val,
                'cosine' : cosine_corr,
                'jaccard' : jaccard_corr
            }
        )
    
    correlation_df = pd.DataFrame(correlation_matrix)
    correlation_df.set_index('drug', inplace=True)
    return correlation_df
