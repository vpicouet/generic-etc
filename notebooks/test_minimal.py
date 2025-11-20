# Test minimal pour identifier le problème de double figure
import matplotlib
matplotlib.use('module://ipympl.backend_nbagg')
import matplotlib.pyplot as plt

# Simuler le contexte notebook
from ipywidgets import Output
from IPython.display import display

out = Output()

print("Test 1: Figure dans Output widget")
with out:
    fig, ax = plt.subplots()
    ax.plot([1,2,3], [1,4,9])
    plt.show()

display(out)
print("Nombre de figures:", len(plt.get_fignums()))
