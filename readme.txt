Project: Ergasia_2_2020-2021

Eleftheria Ellina
A.M. : 1115201800228
Stylianos Psara
A.M. : 1115201800226

1. Ektelesi programmatos

Compile with make and excecute with command below:

FOR SEARCH:
with default values:
./search -i nasd_input.csv -q nasd_query.csv -o outfile1.txt -algorithm LSH -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -o outfile1.txt -algorithm Hypercube -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -o outfile1.txt -algorithm Frechet -metric Discrete -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -o outfile1.txt -algorithm Frechet -metric Continuous -delta 1.0

with values taken from command line:
./search -i nasd_input.csv -q nasd_query.csv -k 4 -L 5 -M 10 -probes 2 -o outfile1.txt -algorithm LSH -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -k 4 -L 5 -M 10 -probes 2 -o outfile1.txt -algorithm Hypercube -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -k 4 -L 5 -M 10 -probes 2 -o outfile1.txt -algorithm Frechet -metric Discrete -delta 1.0
./search -i nasd_input.csv -q nasd_query.csv -k 4 -L 5 -M 10 -probes 2 -o outfile1.txt -algorithm Frechet -metric Continuous -delta 1.0

FOR CLUSTERING:
without the optional printing in output file:
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Classic
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Classic
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Frechet
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment LSH
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Hypercube

with the optional printing in output file:
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Classic -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Classic -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Frechet -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment LSH -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Hypercube -silhouette 1

with the optional printing in output file:
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Classic -complete 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Classic -complete 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Frechet -complete 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment LSH -complete 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Hypercube -complete 1

with the optional printing in output file:
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Classic -complete 1 -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Classic -complete 1 -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Frechet -assignment Frechet -complete 1 -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment LSH -complete 1 -silhouette 1
./cluster -i nasd_input.csv -c cluster.conf -o outfile2.txt -update Mean_Vector -assignment Hypercube -complete 1 -silhouette 1

FOR TEST:
./test

2. Arxeia pou paradidonte

data.hpp
data.cpp

LSH.hpp
LSH.cpp

hypercube.hpp
hypercube.cpp

cluster_main.cc
clustering.hpp
clustering.cpp

Frechet_LSH.cpp
Frechet_clustering.cpp

search_main.cpp
test_main.cpp


data.cpp
	functions: (Apo project 1)
		comp: sigkrinei ta points me metro sigkrisis tin apostasi, einai typou bool
		distance: ipologizei tin eukleidia apostasi meta3i 2 dianismaton kai tin epistrefi
		compare_distance: sigkrini ta points me metro sigkrisis tin apostasi, einai typou bool
		compare_ID_pointer: sigkrini ta points me metro sigkrisis to ID, einai typou bool
	functions: (Apo project 2)
		Discrete_Frechet_Distance: ipologizi tin discrete frechet apostasi me dynamic array
		compare_dist: sigkrini ta points me metro sigkrisis tin apostasi, einai typou bool, pairnei deiktes
		filtering: sigkrini ta simeia tis kampilis me tin dotheisa condition apo tis diafaneies kai diamorfonei tin kainourgia filtrarismeni 
			   kampili analoga
		ContinuousFrechet_distance: dimiourgei ta katalila antikeimena apo tis klaseis ton arxeion tou github pou mas dotheikan kai kalei tin sinartisi
					    distance tou dothentos library

LSH.cpp
	functions: 
		hashPush: vazei ta points se hashtable vasi tis g hashfunction
		hi: ipologizi ti timi h(i) meso tou tipou pou dothike me tis metavlites p, v, t, w
		v_calculator: ipologizi v dianismata meso tis kanonikis katanomis
		inner_product: ipologizei to esoteriko ginomeno v.p
		gp: ipologizi to g bucket meso tou typoy pou mas dothike
		gp_cluster: ipologizi to g bucket gia ta clusters
		modulo:	sinartisi gia modulo pou xrisimopieitai apo tin g gia tin orthi xrisi tou mod
		kNNSearch: pairnei to query ipologizei to g bucket kai sti sinexia sigkrinei ola ta points tou bucket meso ID,
			an iparxei kapio idio tote to apothikevo os ena apo ta kontinotera, stin sinexia gia ta ipolipa kontinotera 
			an iparxoun sigkrini ta points meso eukleidias apostasis, sto telos epistrefi ta N kontinotera.
		range_search: pairnei to query kai ipologizi tis apostasis ton points oi opoies ine entos tou themitou R

hypercube.cpp:
	functions:
		label: ipologismos label gia kathe point pou eisagete sto hashtable
		label_cluster: ipologismos label gia kathe point pou eisagete se cluster
		hashPush: vazei ta points se hashtable vasi ton label
		fi: pairnei timi hi kai thetei bit tixea
		hi: ipologizi ti timi h(i) meso tou tipou pou dothike me tis metavlites p, v, t, w
		v_calculator: ipologizi v dianismata meso tis kanonikis katanomis
		inner_product: ipologizei to esoteriko ginomeno v.p
		hamming_distance: ipologizi to hamming distance meta3i dio label
		kNNsearch: vriski to label gia to query kai vriskei random M points kai epistrefei ta N kontinotera,
			se periptosi pou den vri M simeia allazei bucket kai vriskei ta ipolipa, termatizei otan vri N stoixeia
			i otan ftasei sto max number of buckets.
		range_search: vriski to label gia to query kai vriskei random M points kai epistrefei apo ta M afta pou einai mikrotera apo to R,
			se periptosi pou den vri M simeia allazei bucket kai vriskei ta ipolipa, termatizei otan vri M stoixeia
			i otan ftasei sto max number of buckets. San sxediastiki epilogi exoume vali gia tin efkolia mas sto cluster na perni 
			minR kai maxR, opou stin klisi tou hyper aftonoma perni minR os miden.

clustering.cpp:
	functions:
		init_Push: eisagei ta points sto proto centroid pou thetoume tixea
		Di_first_centroid: arxikopoiei tis apostasis Di gia kathe point apo to proto centroid
		k_means: 3ekinaei me ena centroid sto opoio anatheti ola ta points, kai vriski to epomeno centroid me vasi tou probability typou pou
			mas dothike kai 3anakanei sosto assign kathe point sto katalilo centroid, sinexizi mexri na vri K centroids.
		table_push: anatheti ta points sto sosto centroid otan afta ginoun updated.
		Lloyds: kanei update ta centroids kai anatheti ta points sta sosta cluster, termatizi otan stamatisoun ta centroid na allazun
		new_centroid: vriski to neo centroid meso tou mesou dianismatos
		binary_search: kanei binary search oste na vri tin kontinoteri pithanotita oste na vrethi to epomeno centroid gia tin k-means
		table_Push_reverse: vazi ta points sto sosto centroid gia tin diadikasia tou reverse range search
		reverse_range_search_LSH: pairnei ta clusters meta tin anathesi tu k-means, sti sinexeia vriski me ti range search tis LSH ta simeia
					pou vriskonte entos aktinas R appo to centroid kai ta anathetei, i aktina R diplasiazete kai i epanala4i sinexizi mexri
					na mi ginete alli anathesi sta centroid, otan afto stamatisi anatheti ta ipolipa points sta katalila centorid kai afta pou vriskonte se 2
					i perisotera clusters, anathetonte se 1 katopin erevnas meso tis apostasis gia evresi tou katalilou centroid, sti sinexeia ta centroids ginonte
					updated, kai 3anaginete anathesi ton points sta katallila centroids apo tin arxi, afti i epanali4i stamataei otan stamatisoun na allazoun ta centroids.
		reverse_range_search_HYPER: pairnei ta clusters meta tin anathesi tu k-means, sti sinexeia vriski me ti range search tis hypercube ta simeia
					pou vriskonte entos aktinas R apo to centroid kai ta anathetei, i aktina R diplasiazete kai i epanala4i sinexizi mexri
					na mi ginete alli anathesi sta centroid, otan afto stamatisi anatheti ta ipolipa points sta katalila centorid kai afta pou vriskonte se 2
					i perisotera clusters, anathetonte se 1 katopin erevnas meso tis apostasis gia evresi tou katalilou centroid, sti sinexeia ta centroids ginonte
					updated, kai 3anaginete anathesi ton points sta katallila centroids apo tin arxi, afti i epanali4i stamataei otan stamatisoun na allazoun ta centroids.
		cond: sinartisi pou ipodikniei pia points ine se periostera apo 2 clusters oste na 3eroume na ta periorisoume se 1
		silhouttes: ginete ipologismos ton silhouettes meso tou tipou pou mas dothike
		output: sinartisi pou tiponei sto file
		optional: sinartisi pou tiponi ta optional an zitithoun apo to command line

Frechet_LSH.cpp:
		fucntions:
			DiscreteFrechet_hashPush: kalei tin DiscreteFrechet_create_new_curve kai vazei ta points se hashtable vasi tis g hashfunction
			DiscreteFrechet_create_new_curve: dimiourgei tin grid kampili, diagrafei ta duplicated simeia kai kanei padding gia tin omoiomorfia 
							  olon ton kampilon se d, oste na hasharistei apo tin g tou Frechet_LSH.
			
			DiscreteFrechet_hi: ipologizi ti timi h(i) meso tou tipou pou dothike me tis metavlites p, v, t, w
			DiscreteFrechet_inner_product: ipologizei to esoteriko ginomeno v.p (gia dio distances)
			DiscreteFrechet_gp: ipologizi to g bucket meso tou typoy pou mas dothike
			DiscreteFrechet_gp_cluster: ipologizi to g bucket gia ta clusters
			DiscreteFrechet_kNNsearch: pairnei to query, kai to metatrepi oste na vri ti grid curve tou i opoia tha xrisimopoithi gia na ipologistei to g bucket 
				kai sti sinexia sigkrinei ola ta points tou bucket meso ID,
				an iparxei kapio idio tote to apothikevo os ena apo ta kontinotera, stin sinexia gia ta ipolipa kontinotera 
				an iparxoun sigkrini ta points meso discrete frechet apostasis, sto telos epistrefi ta N kontinotera.
			DiscreteFrechet_range_search: pairnei to query kai ipologizi tis apostasis ton points oi opoies ine entos tou themitou R
			ContinuousFrechet_hashPush: kalei tin ContinuousFrechet_create_new_curve kai vazei ta points se hashtable vasi tis g hashfunction
			ContinuousFrechet_create_new_curve: dimiourgei tin grid kampili, diagrafei ta duplicated simeia kai kanei padding gia tin omoiomorfia 
							  olon ton kampilon se d, oste na hasharistei apo tin g tou Frechet_LSH.
			ContinuousFrechet_hi: ipologizi ti timi h(i) meso tou tipou pou dothike me tis metavlites p, v, t, w
			ContinuousFrechet_inner_product: ipologizei to esoteriko ginomeno v.p (gia dio distances)
			ContinuousFrechet_gp: ipologizi to g bucket meso tou typoy pou mas dothike
			ContinuousFrechet_kNNsearch: pairnei to query, kai to metatrepi oste na vri ti grid curve tou i opoia tha xrisimopoithi gia na ipologistei to g bucket 
				kai sti sinexia sigkrinei ola ta points tou bucket meso ID,
				an iparxei kapio idio tote to apothikevo os ena apo ta kontinotera, stin sinexia gia ta ipolipa kontinotera 
				an iparxoun sigkrini ta points meso continuous frechet apostasis, sto telos epistrefi ta N kontinotera.
		
Frechet_clustering.cpp
		  functions:
			Mean_Curve: i diadikasia ginetai anamesa se dio kampiles, vriskei arxika to dynamic array, sti sinexeia kanei backtracking gia na vrei to mean curve kai sto telos
				    filtrarei tis kampiles oste na exoun megethos d (meiosi tou excecution time kai space).
			Frechet_new_centroid: anathetoume sta filla tou dentrou tis kampiles tou trexon cluster, anevenoume pros ta pano ipologizontas to mean curve mexri na ftasoume sti riza
					     opou einai i teliki kampili pou apotelei to neo centroid.
			Frechet_Di_first_centroid: arxikopoiei tis apostasis Di gia kathe point apo to proto centroid
			Frechet_k_means: 3ekinaei me ena centroid sto opoio anatheti ola ta points, kai vriski to epomeno centroid me vasi tou probability typou pou
				mas dothike kai 3anakanei sosto assign kathe point sto katalilo centroid, sinexizi mexri na vri K centroids.
		 	Frechet_Lloyds: kanei update ta centroids kai anatheti ta points sta sosta cluster, termatizi otan stamatisoun ta centroid na allazun
				new_centroid: vriski to neo centroid meso tou mesou dianismatos
			Frechet_table_Push: anatheti ta points sto sosto centroid otan afta ginoun updated.
			Frechet_table_Push_reverse: vazi ta points sto sosto centroid gia tin diadikasia tou reverse range search
			reverse_range_search_LSH_Frechet: pairnei ta clusters meta tin anathesi tu k-means, sti sinexeia vriski me ti range search tis Frechet_LSH ta simeia
					pou vriskonte entos aktinas R apo to centroid kai ta anathetei, i aktina R diplasiazete kai otan afto stamatisi anatheti ta ipolipa points sta katalila centorid kai afta pou vriskonte se 2
					i perisotera clusters, anathetonte se 1 katopin erevnas meso tis apostasis gia evresi tou katalilou centroid, sti sinexeia ta centroids ginonte
					updated, kai 3anaginete anathesi ton points sta katallila centroids apo tin arxi, afti i epanali4i stamataei otan stamatisoun na allazoun ta centroids.
			cluster_filtering: kanoume filtering stin mean curve mexri to size tis na gini d, i parametros e 3ekinaei me timi 1.0 kai af3anete se kathe iteration 0.1.
			curve_silhouttes: ginete ipologismos ton silhouettes meso tou tipou pou mas dothike

GIA OLES TIS MAIN:
Ginete apla klisi ton sinartiseon kai tipoma sta arxeia.