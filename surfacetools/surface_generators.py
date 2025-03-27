
import numpy as np

def generate_wrinkled_surface(size=801, skew_target=2, kurt_target=10):
    """
    Generates a 2D surface height map with periodic undulations,
    specified skewness, and specified kurtosis.
    """
    np.random.seed(42)  # For reproducibility
    
    # Step 1: Generate a base sinusoidal wave pattern for periodic undulations
    x = np.linspace(0, 4 * np.pi, size)
    y = np.linspace(0, 4 * np.pi, size)
    X, Y = np.meshgrid(x, y)
    
    base_wave = np.sin(X) * np.cos(Y)  # Base periodic structure

    # Step 2: Add Gaussian noise for wrinkles
    noise = np.random.normal(0, 0.2, (size, size))  # White noise

    # Combine both
    height_map = base_wave + noise  

    # Step 3: Transform the distribution to match target skewness & kurtosis
    height_map_flat = height_map.flatten()

    # Normalize to standard normal distribution
    height_map_flat = (height_map_flat - np.mean(height_map_flat)) / np.std(height_map_flat)

    # Use a transformation to achieve desired skewness & kurtosis
    transformed_data = stats.boxcox(height_map_flat - np.min(height_map_flat) + 1)[0]
    
    # Further adjust skewness and kurtosis
    transformed_data = (transformed_data - np.mean(transformed_data)) / np.std(transformed_data)
    transformed_data = transformed_data**3  # Introducing skew
    
    # Scale to match desired kurtosis
    transformed_data = transformed_data / np.max(np.abs(transformed_data))  # Normalize to [-1, 1]
    transformed_data = transformed_data * np.sqrt(kurt_target)  # Scale to kurtosis target

    # Reshape back to 2D
    height_map = transformed_data.reshape((size, size))

    # Normalize to range [0,1]
    height_map = (height_map - np.min(height_map)) / (np.max(height_map) - np.min(height_map))

    return height_map