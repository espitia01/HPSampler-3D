import numpy as np
import matplotlib.pyplot as plt
from HPMove import HPMove
from HPDistance import HPDistance
from HPEnergy import HPEnergy
import time
from HPShow import HPShow
from tqdm import tqdm
import pandas as pd
import os
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
logging.basicConfig(filename='simulation.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

time_dict = {}

def process_protein_length(n):
    print(f"Processing protein length N = {n}")
    local_datasets = []
    local_times = []
    for P in generate_sequences(n):
        start_time = time.time()
        df = HPFold(P)
        end_time = time.time()
        elapsed_time = end_time - start_time
        local_times.append(elapsed_time)
        local_datasets.append(df)
       # Aggregate times for this length
    if n in time_dict:
        time_dict[n].extend(local_times)
    else:
        time_dict[n] = local_times
    return pd.concat(local_datasets)


def translate_to_origin(S):
    return S - S[:, 0].reshape(-1, 1)

def rotate(S, axis, angle):
    rotation_matrix = {
        'x': np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]]),
        'y': np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]]),
        'z': np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    }
    return np.dot(rotation_matrix[axis], S)

def reflect(S, plane):
    reflection_matrix = {
        'xy': np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]),
        'yz': np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        'zx': np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]])
    }
    return np.dot(reflection_matrix[plane], S)

def minimize_form(forms):
    # Sort and find the minimum form lexicographically
    return min(forms, key=lambda x: tuple(x.flatten()))

def canonical_form(S):
    S_translated = translate_to_origin(S)
    
    candidate_forms = []
    
    # Consider rotations around x, y, and z axes
    for axis in ['x', 'y', 'z']:
        for angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
            rotated = rotate(S_translated, axis, angle)
            candidate_forms.append(rotated)
            
            # Also consider the reflected form for each rotation
            for plane in ['xy', 'yz', 'zx']:
                reflected = reflect(rotated, plane)
                candidate_forms.append(reflected)

    # Choose the canonical form as the 'smallest' one
    return minimize_form(candidate_forms)

def hash_conformation(S):
    S_canonical = canonical_form(S)
    return hash(tuple(S_canonical.flatten()))

def generate_sequences(n):
    num_sequences = 2 ** n
    sequences = np.array([[(i >> j) & 1 for j in range(n)] for i in range(num_sequences)], dtype=int)
    return sequences
    
MAX_ITER = 100  # Maximum iterations without progress
iter_without_progress = 0  # Counter for iterations without progress


    
#Initial values
def HPFold(P,Time=200,Temperature=1.5,J=-1.0,Mode=0):
    P = np.array(P)
    N = len(P)
    
    seen_hashes = set()

    #further internal variables
    M = int(np.ceil(np.max([10*4*N,Time])))
    MShow = np.floor(M/10)

    # Initial protein structure
    #S = np.concatenate((np.ones((1,N),dtype-int),np.zeros((1,N),dtype=int)))  #stretched in x
    S = np.zeros((3,N),dtype=int)  # 3D structure
    S[0,1::2] = 1
    S[1,2::2] = 1
    S[:,0] = 0
    D = HPDistance(S)
    # Temperature range
    Temperature = np.max([Temperature,1.5*np.abs(J)])
    Temperature = np.linspace(Temperature,np.abs(J)/20.0,10)

    #show the initial configuration
    SMin = S
    EMin = HPEnergy(P,S,J)
    #HPShow(P,S,EMin,Temperature[0],0)
    #plt.pause(1)
    
    dataset = []

    #the temperature-loop 
    for temp in Temperature:
        S = SMin #start each temperature at the lowest-energy prior structure
        for i in range(M):
            E_i = HPEnergy(P,S,J)
            if N > 2:
                j = np.random.randint(2, N)
            else:
                j = 1
            probability = 0.

            #the Monte-Carlo loop
            #the Monte-Carlo loop
            while probability < np.random.random(1):
                if j < S.shape[1]:
                    s_j = S[:,j]
                else:
                    j = 0
                    s_j = S[:,j]
                Dmin = 0

                #find an allowed move in j without overlap
                while Dmin == 0:
                    s_mv = HPMove(S[:,j], S[:,j-1])
                    S[:,j] = s_mv
                    D_mv = HPDistance(S)
                    #Check for overlap: This will give a 0 value is the D matrix
                    Dmin = np.min(D_mv)
                
                E_new = HPEnergy(P, S, J)
                current_hash = hash_conformation(S)
                if current_hash in seen_hashes:
                    iter_without_progress += 1
                    if iter_without_progress >= MAX_ITER:
                        break
                else:
                    seen_hashes.add(current_hash)
                    iter_without_progress = 0  

                dataset.append((P.copy(), S.copy(), E_new))
                probability = np.exp(-(E_new-E_i)/temp)
                #print(temp, i, j, E_new, E_new-E_i, probability)
                if E_new<EMin:
                    EMin = E_new
                    SMin = S
            
            #a new configuration was found
            #if np.mod(i,MShow)==0: HPShow(P,S,E_new,temp,i)

        #final configuration for that temperature
        #HPShow(P,S,E_new,temp,M)
    

    #the final-final configuration
    df = pd.DataFrame(dataset, columns=['Sequence', 'Structure', 'Energy'])
    min_energy = df['Energy'].min()
    max_energy = df['Energy'].max()
    # Adjust the line where you assign the 'Label' column in your DataFrame as follows:
    df['Label'] = np.where(df['Energy'] == 0.0, 0, np.where(df['Energy'] == min_energy, 1, 0))
    #HPShow(P,SMin,EMin,temp,M)
    return df


if __name__ == "__main__":
    max_length = 50
    all_datasets = []
    time_dict = {}

    logging.info("Starting protein simulations")
    print(f"Starting protein simulations up to length: {max_length}")
    print("Simulation setup complete. Starting now...")

    max_workers = max(1, os.cpu_count() // 2)
    print(f"Max workers: {max_workers}")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for n in range(4, max_length + 1):
            future = executor.submit(process_protein_length, n)
            futures[future] = n
            print(f"Submitting task for protein length: {n}")

        for future in as_completed(futures):
            protein_length = futures[future]
            try:
                result = future.result()
                all_datasets.append(result)
                print(f"Successfully processed protein length: {protein_length}")
                logging.info(f"Completed protein length: {protein_length}")
            except Exception as e:
                logging.error(f"Error processing protein length {protein_length}: {e}")

    combined_df = pd.concat(all_datasets)
    combined_df.to_csv('combined_dataset_parallel.csv', index=False)

    average_time_dict = {n: sum(times) / len(times) for n, times in time_dict.items()}
    logging.info(f"Average Time per Sequence Length: {average_time_dict}")
    print("Average Time per Sequence Length:", average_time_dict)
    print("All protein simulations completed. Saving combined dataset...")
