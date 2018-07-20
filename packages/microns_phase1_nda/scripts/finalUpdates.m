clear tuple
for i=570:579
    tuple.em_id=i;
    insert(nda.EMCell,tuple);
end

%%
clear tuple
tuple.stack_name='2P'
tuple.stack_dimensions=[400 400 310];
tuple.stack_pixels=[512 512 310];
tuple.stack_description='High resolution 2P stack';
insert(nda.Stack,tuple);

tuple.stack_name='Pinky40'
tuple.stack_dimensions=nan;
tuple.stack_pixels=nan;
tuple.stack_description='Pinky40 EM Stack';
insert(nda.Stack,tuple);

%%
for i=1:length(c)
    clear tuple
    tuple.em_id=c(i).id;
    tuple.stack_name='2P';
    tuple.soma_x=c(i).xyz(1);
    tuple.soma_y=c(i).xyz(2);
    tuple.soma_z=c(i).xyz(3);
    insert(nda.Soma,tuple);
    
    if c(i).inEM
        tuple.em_id=c(i).id;
        tuple.stack_name='Pinky40';
        tuple.soma_x=c(i).xyzEM(1);
        tuple.soma_y=c(i).xyzEM(2);
        tuple.soma_z=c(i).xyzEM(3);
        insert(nda.Soma,tuple);
    end
end
    
    
    