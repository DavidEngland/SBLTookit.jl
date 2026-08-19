# src/Catalog.jl
module Catalog

using JSON3

struct DatasetMeta
    id::String
    name::String
    status::String
    priority::String
    ingestion_source::String
    primary_papers::Vector{String}
end

function load_catalog(json_path::String)::Dict{String, DatasetMeta}
    raw = read(json_path, String)
    data = JSON3.read(raw)
    catalog = Dict{String, DatasetMeta}()

    for ds in data.datasets
        id = String(ds.id)
        meta = DatasetMeta(
            id,
            String(ds.name),
            String(ds.status),
            String(ds.priority),
            String(ds.data_assets.ingestion_source),
            String.[p for p in ds.literature_context.primary_papers]
        )
        catalog[id] = meta
    end
    return catalog
end

end # module