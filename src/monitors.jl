# Function to parse column headers
function get_cols(head::Vector{String})
    cols = String[]
    
    for i in 1:length(head)
        io = 1
        for j in 1:length(split(head[1]))
            if i == 1
                if j == 1
                    push!(cols, strip(head[i][1:12]))
                    io = 13
                else
                    push!(cols, strip(head[i][io:io + 13]))
                    io += 14
                end
            else
                if j == 1
                    w = strip(head[i][1:12])
                    if w != ""
                        cols[j] *= "_" * w
                    end
                    io = 13
                else
                    w = split(strip(head[i][io:io + 13]))
                    for k in 1:length(w)
                        if w[k] != ""
                            cols[j] *= "_" * w[k]
                        end
                    end
                    io += 14
                end
            end
        end
    end

    # Trim whitespace from each column name
    return [strip(v) for v in cols]
end

# Function to parse data
function get_data(dat::Vector{String}, ns::Int, nv::Int)
    dat_np = zeros(Float64, ns, nv)
    for i in 1:ns
        step = split(dat[i])
        for j in 1:nv
            dat_np[i, j] = parse(Float64, step[j])
        end
    end
    return dat_np
end

# Main function to read and convert to DataFrame
function mon2df(fn::String)
    lines = readlines(fn)
    
    # Determine the number of timesteps and header length
    ns = parse(Int, split(lines[end])[1])
    nh = length(lines) - ns
    
    # Split into header and data
    head = lines[1:nh-1]
    dat = lines[nh:end]
    
    # Get columns and data
    cols = get_cols(head)
    nv = length(cols)
    dat_np = get_data(dat, ns, nv)
    
    # Create a DataFrame
    df = DataFrame(dat_np, Symbol.(cols))
    return df
end


function mean_noinf(arr)
    m=0.0
    c=0
    for i in eachindex(arr)
        if abs(arr[i]) != Inf
            m+=arr[i]
            c+=1
        end
    end
    return m/c
end
