import numpy as np
import matplotlib.pyplot as plt
from HPMove import HPMove
from HPDistance import HPDistance
from HPEnergy import HPEnergy
import time
from HPShow import HPShow
import pandas as pd

time_dict = {}


def translate_to_origin(S):
    return S - S[:, 0].reshape(-1, 1)

def rotate(S, angle):
    rotation_matrix = np.array([[np.cos(angle), -np.sin(angle)],
                                [np.sin(angle), np.cos(angle)]])
    return rotation_matrix.dot(S)

def reflect(S):
    # Reflect over the y-axis
    return np.array([[S[0, i], -S[1, i]] for i in range(S.shape[1])]).T

def minimize_form(forms):
    # Sort and find the minimum form lexicographically
    return min(forms, key=lambda x: tuple(x.flatten()))

def canonical_form(S):
    S_translated = translate_to_origin(S)
    
    candidate_forms = []
    
    # Consider 0, 90, 180, 270 degree rotations
    for angle in [0, np.pi/2, np.pi, 3*np.pi/2]:
        rotated = rotate(S_translated, angle)
        candidate_forms.append(rotated)
        
        # Also consider the reflected form for each rotation
        reflected = reflect(rotated)
        candidate_forms.append(reflected)

    # Choose the canonical form as the 'smallest' one
    return minimize_form(candidate_forms)

def hash_conformation(S):
    S_canonical = canonical_form(S)
    return hash(tuple(S_canonical.flatten()))

def generate_sequences(n):
    if n == 1:
        return [[0], [1]]
    else:
        sequences = []
        for sequence in generate_sequences(n - 1):
            sequences.append(sequence + [0])
            sequences.append(sequence + [1])
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
    S = np.zeros((2,N),dtype=int)  # staircase
    S[0,1::2] = 1
    S[1,2::2] = 1
    S[:,0] = 0
    
    # Temperature range
    Temperature = np.max([Temperature,1.5*np.abs(J)])
    Temperature = np.linspace(Temperature,np.abs(J)/20.0,10)

    #show the initial configuration
    SMin = S
    EMin = HPEnergy(P,S,J)
    HPShow(P,S,EMin,Temperature[0],0)
    plt.pause(1)
    
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
                
                E_new = HPEnergy(P,S,J)
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
            if np.mod(i,MShow)==0: HPShow(P,S,E_new,temp,i)

        #final configuration for that temperature
        HPShow(P,S,E_new,temp,M)
    

    #the final-final configuration
    df = pd.DataFrame(dataset, columns=['Sequence', 'Structure', 'Energy'])
    min_energy = df['Energy'].min()
    max_energy = df['Energy'].max()
    df['Label'] = np.where(df['Energy'] == min_energy, 1, (np.where(df['Energy'] != min_energy, 0, np.nan)))
    HPShow(P,SMin,EMin,temp,M)
    return df

max_length = 4

all_datasets = []

for n in range(4, max_length + 1):
    for P in generate_sequences(n):
        start_time = time.time()  # Record start time
        df = HPFold(P)
        end_time = time.time()  # Record end time
        elapsed_time = end_time - start_time  # Calculate elapsed time

        # Store elapsed time in the dictionary
        if n in time_dict:
            time_dict[n].append(elapsed_time)
        else:
            time_dict[n] = [elapsed_time]
        all_datasets.append(df)

combined_df = pd.concat(all_datasets)
combined_df.to_csv('combined_dataset_2.csv', index=False)

average_time_dict = {n: sum(times)/len(times) for n, times in time_dict.items()}
print("Average Time per Sequence Length:", average_time_dict)
