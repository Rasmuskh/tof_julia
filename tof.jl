using DataFrames
using Plots

function basic_framer(filepath, threshold)
    #Get the number of lines
    function get_lines(path)
        open(path) do file
            nlines=0
            for ln in eachline(file)
                nlines += 1
            end
        end
        return nlines
    end

    nlines = get_lines("$filepath")
    nevents = trunc(Int,nlines/8)
    #get the timestamp and the samples
    samples = Vector(undef, trunc(Int, nevents))
    timestamp = Vector(undef, trunc(Int,nevents))
    line_index = 0
    event_index = 1
    open("$filepath") do file
        for line in eachline(file)
            line_index += 1
            if line_index%8 == 6
                print(trunc(Int, line_index/nlines*100), "%\r")
                #specify type
                timestamp[event_index]=parse(Int, split(line, ": ")[2])
            elseif line_index%8==0
                #specify type
                samples[event_index] = [parse(Int, w) for w in split(strip(line), " ")]
                samples[event_index] .-= trunc(Int, (sum(samples[event_index][1:20])/20))
                if minimum(samples[event_index])>threshold
                    event_index+=1
                elseif minimum(samples[event_index])<-threshold
                    samples[event_index] .*= -1
                    event_index+=1
                end
            end
        end
    end
    samples = samples[1:event_index-1]
    timestamp = timestamp[1:event_index-1]
    return DataFrame([samples, timestamp], [:samples, :timestamp])
end


function cfd(samples, frac)
    peak = maximum(samples)
    crossing = 0
    for i in 1:length(samples)
        if samples[i] >= frac*peak
            crossing = i
            break
        end
    end
    return crossing
end


function advanced_framer(frame, cfd_frac=0.5)
    nTimeResets = 0
    timestamp = Vector(undef, size(frame, 1))
    refpoint = Vector(undef, size(frame, 1))
    for n in 1:size(frame, 1)
        if n%10 == 0
            print(trunc(Int, n/size(frame,1)*100), "%\r")
        end
        if n > 1
            if frame.timestamp[n]<frame.timestamp[n-1]
                nTimeResets += 1
            end
        end
        timestamp[n] = (frame.timestamp[n]+nTimeResets*2147483647)
        refpoint[n] = cfd(frame.samples[n], cfd_frac)
    end

    return DataFrame([timestamp, refpoint], [:timestamp, :refpoint])
end

function tof_spectrum(neutron, gamma, fac=16)
    tolerance = 100
    ymin = 0
    Thist = zeros(Int16, 200, 1)
    for ne in 1:size(neutron, 1)
        if ne%10==0
            print(trunc(Int,ne/size(neutron,1)*100), "%\r")
        end
        for y in 1:size(gamma, 1)
            Delta = (fac*neutron.timestamp[ne]+neutron.refpoint[ne]) - (fac*gamma.timestamp[y]+gamma.refpoint[y])
            if Delta>tolerance
                ymin=y
            end
            if -tolerance < Delta < tolerance
                Thist[tolerance+trunc(Int16, Delta)]+=1
            elseif Delta<-tolerance
                break
            end
        end
    end
    return Thist
end



@time begin
    #get neutron frame
    print("generating basic neutron frame\n")
    #basic_neutron = basic_framer("data/2018-07-27/4minutes2018-07-27-N10-G10ch1.txt", 10)
    basic_neutron = basic_framer("testneutron.txt", 20)
    print("generating advanced neutron frame\n")
    adv_neutron = advanced_framer(basic_neutron)
    basic_neutron = nothing
    GC.gc()
    #get gamma frame
    print("generating basic gamma frame\n")
    #basic_gamma = basic_framer("data/2018-07-27/4minutes2018-07-27-N10-G10ch1.txt", 10)
    basic_gamma = basic_framer("testgamma.txt", 20)
    print("generating advanced gamma frame\n")
    adv_gamma = advanced_framer(basic_gamma)
    basic_gamma = nothing
    GC.gc()
    #get tof spectrum
    print("generating tof spectrum\n")
    spectrum = tof_spectrum(adv_neutron, adv_gamma)
    #plot stuff
    unicodeplots()
    plot(spectrum)
    plotly()
    plot(spectrum)
end

