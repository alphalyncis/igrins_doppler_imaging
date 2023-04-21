import matplotlib.pyplot as plt
import numpy as np
import paths

print("import ok.")
# Generate some data
random_numbers = np.random.randn(100, 10)

# Plot and save
fig = plt.figure(figsize=(7, 6))
plt.plot(random_numbers)
plt.xlabel("x")
plt.ylabel("y")
fig.savefig(paths.figures / "random_numbers.pdf", bbox_inches="tight", dpi=300)
print("Figure saved.")