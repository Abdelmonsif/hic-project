#python graph_preprocess.py -data_dir ../../test_data/test_0.gexf -edge_dir ../../test_data/test_0.h5 -node_dir ../../test_data/test_0.csv
#python graph_preprocess.py -data_dir ../../test_data/test_1.gexf -edge_dir ../../test_data/test_1.h5 -node_dir ../../test_data/test_1.csv
#python graph_preprocess.py -data_dir ../../test_data/test_2.gexf -edge_dir ../../test_data/test_2.h5 -node_dir ../../test_data/test_2.csv
#python graph_preprocess.py -data_dir ../../test_data/test_3.gexf -edge_dir ../../test_data/test_3.h5 -node_dir ../../test_data/test_3.csv
#python graph_preprocess.py -data_dir ../../test_data/test_4.gexf -edge_dir ../../test_data/test_4.h5 -node_dir ../../test_data/test_4.csv

#python recover_gexf.py -edge_dir ../../test_data/test_0.h5 -node_dir ../../test_data/test_0.csv -recover_dir ../../test_data/test_0_recovered.gexf
#python recover_gexf.py -edge_dir ../../test_data/test_1.h5 -node_dir ../../test_data/test_1.csv -recover_dir ../../test_data/test_1_recovered.gexf
#python recover_gexf.py -edge_dir ../../test_data/test_2.h5 -node_dir ../../test_data/test_2.csv -recover_dir ../../test_data/test_2_recovered.gexf
#python recover_gexf.py -edge_dir ../../test_data/test_3.h5 -node_dir ../../test_data/test_3.csv -recover_dir ../../test_data/test_3_recovered.gexf
#python recover_gexf.py -edge_dir ../../test_data/test_4.h5 -node_dir ../../test_data/test_4.csv -recover_dir ../../test_data/test_4_recovered.gexf

#python graph_isomorphism_test.py -data_old ../../test_data/test_0.gexf -data_new ../../test_data/test_0_recovered.gexf
#python graph_isomorphism_test.py -data_old ../../test_data/test_1.gexf -data_new ../../test_data/test_1_recovered.gexf
#python graph_isomorphism_test.py -data_old ../../test_data/test_2.gexf -data_new ../../test_data/test_2_recovered.gexf
#python graph_isomorphism_test.py -data_old ../../test_data/test_3.gexf -data_new ../../test_data/test_3_recovered.gexf
python graph_isomorphism_test.py -data_old ../../test_data/test_4.gexf -data_new ../../test_data/test_4_recovered.gexf