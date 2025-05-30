echo "Running for iron"
#./RunEventLoop_plasticSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 26 minervame1A > plastic_iron_wgt1.txt
#./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 26 minervame1A > iron_wgt1.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 26 minervame1A > 1_26_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 2 26 minervame1A > 2_26_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 3 26 minervame1A > 3_26_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 5 26 minervame1A > 5_26_phys.txt
echo "Running for lead"
#./RunEventLoop_plasticSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 82 minervame1A > plastic_lead_wgt1.txt
#./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 82 minervame1A > lead_wgt1.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 1 82 minervame1A > 1_82_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 2 82 minervame1A > 2_82_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 3 82 minervame1A > 3_82_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 4 82 minervame1A > 4_82_phys.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 5 82 minervame1A > 5_82_phys.txt
echo "Running for carbon"
#./RunEventLoop_plasticSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 3 6 minervame1A > plastic_carbon_wgt1.txt
#./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 3 6 minervame1A > carbon_wgt1.txt
./RunEventLoop_physSB /minerva/data/users/fakbar/NukeHists/v1_NuMAD/ 3 6 minervame1A > 3_6_phys.txt
echo "Done Running"
