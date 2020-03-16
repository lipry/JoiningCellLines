from src.config.config import config
from src.tasks.intersection_mode import intersection_mode_exec
from src.tasks.preprocess_mode import preprocess_mode_exec

exp = config['mode']

if exp == "inter":
    c = config['intersection_mode']
    intersection_mode_exec(c)

if exp == "preprocess":
    c = config['preprocess']
    preprocess_mode_exec()
