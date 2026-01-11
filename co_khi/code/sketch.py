import numpy as np
import matplotlib.pyplot as plt

# Tạo dữ liệu x
x = np.linspace(0.02, 0.15, 400)

# Định nghĩa hàm f(x)
y = 49 / (75 * x)

# Vẽ đồ thị
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Đồ thị f(x) = 49 / (75x)')
plt.grid(True)

plt.show()
