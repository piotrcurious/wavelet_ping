import numpy as np
import matplotlib.pyplot as plt

def plot_signal(filename, title, output_png):
    data = np.fromfile(filename, dtype=np.float32)
    plt.figure(figsize=(10, 4))
    plt.plot(data)
    plt.title(title)
    plt.xlabel("Samples")
    plt.ylabel("Amplitude")
    plt.grid(True)
    plt.savefig(output_png)
    plt.close()
    print(f"Saved {output_png}")

if __name__ == "__main__":
    plot_signal("recorded_buffer.bin", "Recorded Buffer (with noise and echoes)", "recorded_buffer.png")
    plot_signal("initial_xcorr.bin", "Initial Cross-Correlation", "initial_xcorr.png")
    plot_signal("final_buffer.bin", "Final Buffer (after echo subtraction)", "final_buffer.png")
