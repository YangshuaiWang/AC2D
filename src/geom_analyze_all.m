function geom = geom_analyze_all(geom)
    geom = geom_analyze(geom, 1);
    geom = gqc23_prep_geom(geom);
    geom = geom_2dtri_h(geom);
end