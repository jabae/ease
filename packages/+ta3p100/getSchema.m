function obj = getSchema
persistent schemaObject
if isempty(schemaObject)
    schemaObject = dj.Schema(dj.conn, 'ta3p100', 'microns_ta3p100');
end
obj = schemaObject;
end
