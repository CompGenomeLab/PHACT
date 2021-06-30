#Imports from Jinja2
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
import os
import sys

#Import YAML from PyYAML
import yaml

#Load data from YAML file into Python dictionary
cluster_name = sys.argv[1]
if (cluster_name == "sabanci"):
	config = yaml.load(open('config_sabanci.yml'))
else:
	config = yaml.load(open('config_truba.yml'))

#Load Jinja2 template
env = Environment(loader = FileSystemLoader('./'), trim_blocks=True, lstrip_blocks=True)
template = env.get_template('config_template.yml')

# get current directory to set workdir
root_path=os.getcwd()
os.chdir('../../../')
pwd=os.getcwd()
template.globals['workdir'] = pwd
os.chdir(root_path)

# determine the proteins to be run
with open(config['cluster']['query_file']) as file_in:
    lines = []
    for line in file_in:
        lines.append(line.strip())
template.globals['query_ids'] = lines

# set archived path for check
template.globals['archived_path'] = config['cluster']['archived_path']


#Render template using data and print the output
print(template.render(config))
