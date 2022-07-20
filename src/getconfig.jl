using Pkg.TOML

getconfig(key) = getconfig()[key]

function getconfig()
    configfile = joinpath(homedir(), ".GlobalEnergyGIS_config")
    if !isfile(configfile)
        error("Configuration file missing, please run saveconfig(datafolder, uid, api_key) first. See GlobalEnergyGIS README.")
    end
    return TOML.parsefile(configfile)
end
