function TlineFDTD(C,L,d,nCells,cfln,simtime,Vg,Rg,RL,vpts,ipts)
#=
This function models a 1D transmission line using the FDTD algorithm. The
transmission line is excited via a Thevenin source and feeding a resistive
load.
--------------------------------------------------------------------------
Inputs:
C,L = Capacitance and Inductance per unit length
d = length of the transmission line
ncells = number of desired cells
cfln = the CFL number
simtime = total simulation time
Vg = source voltage function Vg(t)
Rg,RL = source and load resistances
vpts,ipts = list of line voltage and current sample points for output
=#

    # Calculate Line Impedance and Travel Speed
    Zo = sqrt.(L./C)
    vo = 1.0./sqrt.(L.*C)

    # Calculate discretization parameters
    dx = sum(d)./nCells
    dt = cfln .* dx ./ maximum(vo)
    nx = nCells + length(RL)
    nt = Int(floor(simtime./dt))

    # Zero initial conditions
    V = zeros(nx,1)
    I = zeros(nx-1,1)

    # Compute Source Parameters
    b1 = C[1].*dx.*0.5./dt
    b2S = 0.5./Rg
    c1S = 1.0./(b1+b2S)
    c2S = b1-b2S
    # Compute Load Parameters
    if typeof(RL) == Array{Int64,1}
        RL = convert(Array{Float64},RL)
    end
    if length(RL) > 1
        for ind = 1:length(RL)-1
            if RL[ind] != 0
                parallel = 1/RL[ind]+1/Zo[ind+1]
                RL[ind] = 1/parallel
            end
        end
    end
    b2L = 0.5./RL
    c1L = 1.0./(b1+b2L)
    c2L = b1-b2L
    # Compute Line Parameters
    cv = dt./(C.*dx)
    ci = dt./(L.*dx)

    # Initialize Probe Vectors
    iProbe = zeros(nt,length(ipts))
    vProbe = zeros(nt,length(vpts))

    # Load Locations
    Rloc = Int.(d./dx + 1)
    if length(Rloc)>1
        for i = 2:length(Rloc)
            Rloc[i] = Rloc[i]+Rloc[i-1]
        end
    end


    # Time Loop
    for n = 1:nt
        # Update voltages
        #= vupdate     # update interior voltages =#
        for k=2:Rloc[1]-1
            V[k] = V[k] - cv[1].*(I[k]-I[k-1])
        end
        if length(Rloc) > 1
            for locInd = 1:(length(Rloc)-1)
                for k=Rloc[locInd]+1:Rloc[locInd+1]-1
                    V[k] = V[k] - cv[locInd].*(I[k]-I[k-1])
                end
            end
        end

        #= vsupdate    # update source voltage =#
        Vs = Vg(n,dt)
        if Rg>0
            V[1] = c1S*(c2S*V[1]-I[1]+Vs/Rg) 
        else
            V[1] = Vs
        end

        #= vLupdate    # update load voltage =#
        if length(Rloc) > 1
            for locInd = 1:length(Rloc)
                if RL[locInd] == 0
                    V[Rloc[locInd]] = 0.
                else
                    V[Rloc[locInd]] = c1L[locInd]*(c2L[locInd]*V[Rloc[locInd]]+I[Rloc[locInd]-1])
                end
            end
        else
            if RL[end]==0
                V[nx] = 0.
            else
                V[nx] = c1L*(c2L*V[nx]+I[nx-1])
            end
        end

        # Update currents
        for k=1:Rloc[1]-1
            I[k] = I[k] - ci[1]*(V[k+1]-V[k])
        end
        if length(Rloc) > 1
            for locInd = 1:(length(Rloc)-1)
                for k=Rloc[locInd]:Rloc[locInd+1]-1
                    I[k] = I[k] - ci[locInd]*(V[k+1]-V[k])
                end
            end
        end

        # Create Output
        for i = 1:length(ipts)
            iProbe[n,i] = I[ipts[i]]
        end
        for v = 1:length(vpts)
            vProbe[n,v] = V[vpts[v]]
        end
    end
    return vProbe, iProbe
end
