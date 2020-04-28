#############################################
########### Output Functions ################
#############################################
function generate_output_data(CRS, savedir)
    unique_name = savedir*"/"*string(hash(CRS.parameters))*Random.randstring()*".csv"
    parameter_df = DataFrame()
    for key in keys(CRS.parameters)
        parameter_df[key] = CRS.parameters[key]
    end
    parameter_df[:save_name] = unique_name
    return unique_name, parameter_df 
end


function InitializeOutput(CRS)
    outDF = DataFrame()
    ### Unclear if this always alines correctly. Could use join(..., by= ) to ensure that
    outDF[Symbol("mID")] = 1:length(CRS.molecule_list)
    outDF[Symbol("molecule")] = CRS.molecule_list
    return outDF
end


function df2json(df::DataFrame)
  len = length(df[:,1])
  indices = names(df)
  jsonarray = [Dict([string(index) => df[index][i] for index in indices])
               for i in 1:len]
  return JSON.json(jsonarray)
end

function writejson(path::String,df::DataFrame)
  open(path,"w") do f
    write(f,df2json(df))
  end
end
