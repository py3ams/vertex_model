function unit_vector = calculate_unit_vector(node_1_position,node_2_position)

adj = node_2_position(1) - node_1_position(1);
opp = node_2_position(2) - node_1_position(2);

hyp = hypot(adj,opp);

unit_vector = 1/hyp*[adj opp];