import warnings
from utils_GREEDY_vertices_SUBSTITUTE import *
warnings.simplefilter('ignore')


case = int(sys.argv[1]) # {1,2}
dim = int(sys.argv[2]) # {1,2}
fe = int(sys.argv[3]) # {0,1}

# Alternatives
ranking_list = list(permutations((0, 1, 2, 3)))
index = 0
pref= ranking_list[index]

import os
_data = sys.argv[4] if len(sys.argv) > 4 else os.environ.get(
    'FISH_DATA', '../../data/fish/fish_J4_K2.csv')
_df = pd.read_csv(_data)
print('characteristics from', _data, flush=True)
print(_df.to_string(index=False), flush=True)
features = _df.drop(columns=['alt_id']).to_numpy(dtype=float)

_SEED = int(os.environ.get('SEED_BASE', 20240808)) + 100 * dim + 10 * fe + case
np.random.seed(_SEED)
print('seed', _SEED, flush=True)
  
if dim == 2:
    features = np.concatenate((features,np.square(features)),axis=1)

if fe == 1:
    fixed_effects_dict = [[x,y,z] for x in range(-10,11) for y in range(-10,11) for z in range(-10,11) ]
    fixed_effects= np.array(fixed_effects_dict[index])
    fixed_effects = np.append(fixed_effects,0)
    fixed_effects = np.array(fixed_effects)
else:
    fixed_effects = np.zeros(features.shape[0])

if case==1:
    x=1
    x_r=0
elif case==2:
    x=2
    x_r=0
elif case==3:
    x=3
    x_r=0
elif case==4:
    x=0
    x_r=1
elif case==5:
    x=2
    x_r=1
elif case==6:
    x=3
    x_r=1
elif case==7:
    x=0
    x_r=2
elif case==8:
    x=1
    x_r=2
elif case==9:
    x=3
    x_r=2
elif case==10:
    x=0
    x_r=3
elif case==11:
    x=1
    x_r=3
elif case==12:
    x=2
    x_r=3
####################################
########Greedy Algorithm############
####################################

greedy_result = process_greedy(pref,features,fixed_effects,x,x_r)
d = {"d":dim,"FE":fe,"chosen":x,"removed":x_r}
print(d)
print(greedy_result)

_out = os.environ.get('OUT_DIR', 'output/raw')
os.makedirs(_out, exist_ok=True)
file_name = _out + "/subs_d" + str(dim) + "_fe" + str(fe) + "_" + str(x) + "_" + str(x_r) + ".txt"
np.savetxt(file_name, greedy_result[1][-1:])
