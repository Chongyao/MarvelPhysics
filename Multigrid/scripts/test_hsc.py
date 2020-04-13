import json
import os
from datetime import datetime

def V_cycle(p):
    p['pre_sm'] = True
    p['post_sm'] = True
    p['diag_PD'] = False
    p['num_mu'] = 1
    p['num_V'] = 30
    p['gs_itrs'] = 4

def W_cycle(p):
    p['pre_sm'] = True
    p['post_sm'] = True
    p['diag_PD'] = False
    p['num_mu'] = 2
    p['num_V'] = 30
    p['gs_itrs'] = 4

def HBF(p):
    p['pre_sm'] = False
    p['post_sm'] = False
    p['diag_PD'] = True
    p['num_mu'] = 1
    p['num_V'] = 1
    p['gs_itrs'] = 10
    
def Hybrid(p):
    p['pre_sm'] = True
    p['post_sm'] = True
    p['diag_PD'] = True
    p['num_mu'] = 1
    p['num_V'] = 1
    p['gs_itrs'] = 1
    

DT = datetime.now().strftime("%b-%d-%Y-%H-%M-%S");
print('date time is ' + DT)

paras = {}

common = {}

exe_path = '../../build/bin/'
data_path = '../data/'
result_path = '../results/'

exe = exe_path + 'test_hsc'
model_name = 'cube_hex'

common['mesh'] = os.path.join(data_path, model_name,  model_name + '.vtk')
common['type'] = 'hex'
common['top_fixed'] = os.path.join(data_path, model_name, model_name + '_top.csv')
common['bottom_fixed'] = os.path.join(data_path, model_name,  model_name + '_bottom.csv') 
common['outdir'] = os.path.join(result_path, 'hsc',model_name, DT)
paras['common'] = common

physics = {}
physics['k'] = 10
physics['w_pos'] = 1e6
paras['physics'] = physics

paras['coarest_num'] = 10
HBF(paras)
paras['PD'] = False

os.makedirs(common['outdir'])
jsonfile = common['outdir'] + '/config.json'
with open(jsonfile, 'w') as f:
    json.dump(paras, f)

cmd = exe + ' ' +  jsonfile
print(cmd)
os.system(cmd)
