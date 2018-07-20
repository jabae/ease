function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'nda', 'microns_nda');
end
obj = schemaObject;
end
