# System Identification using Kernel-Based Regularization

This repository contains a MATLAB implementation of system identification techniques using kernel-based methods. The project focuses on estimating the impulse response of a nonparametric Output Error (OE) model by employing various kernel structures and analyzing the impact of hyperparameters on estimation accuracy.

---

## **Project Overview**

The objective is to identify a system characterized by a nonparametric OE model:


$$y(t) = F_{\theta_0}(z)u(t-1) + e_0(t)$$

where $e_0$ is White Gaussian Noise (WGN) with a standard deviation $\sigma_0 = 0.58$. The system is approximated using an ARX model structure with a length of $n_B = 40$.

The estimation utilizes a Bayesian framework where the parameter vector $\theta$ is assumed to follow a Gaussian distribution $\mathcal{N}(0, K)$, where $K$ is a kernel-based covariance matrix.

### **Key Features**

* **Kernel Implementations**: Supports **Tuned/Correlated (TC)**, **Diagonal (DI)**, and **Stable Spline (SS)** kernel structures.


* **Hyperparameter Analysis**: Evaluates performance across different values of $\lambda$ (scaling factor) and $\beta$ (decay rate).


* **Regularized Estimation**: Implements the `arx` estimation method with regularization matrices derived from the inverse of the kernels.


* **Confidence Intervals**: Calculates and plots 95% confidence intervals to provide a "certificate" of estimation performance.



---

## **Mathematical Framework**

### **Kernel Structures**

Three distinct kernel matrices $K \in \mathbb{R}^{n_B \times n_B}$ are implemented:

* 
**TC (Tuned/Correlated)**: $K(i,j) = \lambda \beta^{\max(i,j)}$.


* 
**DI (Diagonal)**: A diagonal matrix where $K(i,i) = \lambda \beta^i$.


* 
**SS (Stable Spline)**: A correlation structure designed to model smooth, fading impulse responses.



### **Estimation & Confidence Intervals**

The kernel-based estimate $\hat{\theta}_K$ is computed using the dataset $(y^N, u^N)$. To validate the model, the posterior covariance matrix $\hat{P}$ is calculated as:


$$\hat{P} = K - K\Phi^T(\Phi K\Phi^T + \hat{\sigma}_K^2 I)^{-1}\Phi K$$

The 95% confidence interval for each impulse response coefficient is defined by $1.96\sqrt{(\hat{P})_{k+1,k+1}}$.

---

## **Repository Content**

* `Lab304.m`: The primary MATLAB script. It handles data loading, Toeplitz matrix construction ($\Phi$), kernel generation, and the iterative testing of hyperparameters.


* `data.mat`: Contains the input signal $u$, output signal $y$, and the true impulse response coefficients $\theta_0$.



---

## **Usage**

1. Ensure the **System Identification Toolbox** is installed in MATLAB.
2. Place `data.mat` in the same directory as the script.
3. Run `Lab304.m`.


---

## **Conclusions**

The implementation highlights several critical aspects of regularized system identification:

* **Hyperparameter Sensitivity**: Large $\lambda$ values significantly increase the variance of the estimate if the noise level is high.


* **Kernel Selection**: Different kernel types (TC, DI, SS) impose different prior beliefs on the system's smoothness and decay, affecting the bias-variance tradeoff.


* **Reliability**: The inclusion of confidence intervals allows for a visual assessment of where the model is most certain about the system dynamics.
