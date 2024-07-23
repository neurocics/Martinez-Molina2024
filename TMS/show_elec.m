function show_elec(GLAN)

if isfield(GLAN,'chanlocs')
    chanlocs = GLAN.chanlocs;
else
    chanlocs = GLAN;
end

figure_lan
topoplot_lan([],GLAN.chanlocs,'style','blank','electrodes','labelpoint');  

