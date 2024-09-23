# Euler_exact_Riemann matlab code

## Highlight

Compute the exact Riemann solution of the 1D compressible Euler equations (only for ideal polytopic gas) and store key physical quantities so that the exact solution at any point can be retracted using these quantities. 
Good and easy-to-use enough for plotting, code verifying and teaching. 
Support vacuum and zero-temperature/pressure inputs and outputs. Available to output the free boundary speed in the presence of vacuum. Quite rubost. 

## Structure

1. `Euler_exact_Riemann_core.m`: The core function to compute some parameters that cahracterizes the solution, which can then used to analytically compute the solution at any point.
2. `Euler_exact_Riemann_sample.m`: The sampling function, which evaluates the exact solution at user-input positions (mainly for the purpose of plotting).
3. `test.m`: A simple script containing several different test cases.
4. Folder `examples`: Contain some figure outputs that should be expected from running the `test.m`. 

## Examples

![Lax shock tube](https://github.com/user-attachments/assets/ec867b7c-7e61-4423-a928-6ce522fda0d1)
![Sod shock tube](https://github.com/user-attachments/assets/88cad20f-214d-46f6-a18b-272ded49273d)
![LeBlanc shock tube](https://github.com/user-attachments/assets/671029d0-2d4a-459a-ae57-4e8514d9ab26)
![Double rarefaction I](https://github.com/user-attachments/assets/ce947532-983e-4478-adb3-bae981eabb76)
![Double shock](https://github.com/user-attachments/assets/2a92031b-af13-4fa7-9d97-6b7cd1a08eb6)
![Double rarefaction II](https://github.com/user-attachments/assets/8b607238-d0a6-4eec-a50c-a187c05e7530)
![Vacuum middle state](https://github.com/user-attachments/assets/91527b59-34fd-4706-a68a-30c995157f46)
![Vacuum left state](https://github.com/user-attachments/assets/c9c86cf1-27d0-4064-ae4a-e101ff164b16)
