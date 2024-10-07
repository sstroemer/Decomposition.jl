Base.show(io::IO, model::Benders.DecomposedModel) = print(io, "DecomposedModel [algorithm=Benders]: $(model.name)")
