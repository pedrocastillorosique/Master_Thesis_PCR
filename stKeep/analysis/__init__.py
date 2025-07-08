from __main__ import *

import sys
import importlib

module_name_1 = 'analysis.Preprocess_Cell_module'
module_name_2 = 'analysis.Cell_model'

for module_name in [module_name_1, module_name_2]:
    if module_name in sys.modules:
        del sys.modules[module_name]  # Remove from cache
    importlib.import_module(module_name)  # Fresh import