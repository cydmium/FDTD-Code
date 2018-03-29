function Vg(n,dt)
#= 
This function defines the source voltage for the 1D FDTD excitation. It
        represents a trapezoidal pulse with 40 ps rise and fall time with a 100
        ps duration
-------------------------------------------------------------------------------
Inputs:
n = discrete index of time
dt = amount of time between each index
=#

    # Calculate Time Value
    t = n*dt

    # Peak of trapezoid
    if t>=2.6e-10 && t<=3.6e-10
        return 1
    # Rise of trapezoid
    elseif t<2.6e-10 && t>2.2e-10
        return 1*(t-0.22e-9)/0.04e-9
    # Fall of trapezoid
    elseif t<4e-10 && t>3.6e-10
        return 1-1*(t-0.36e-9)/0.04e-9
    # Otherwise 0
    else
        return 0
    end
end
