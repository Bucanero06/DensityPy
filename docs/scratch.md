https://chat.openai.com/share/3da525a7-81f8-46b8-9361-d55ba2151501

sim/ -> PULSE not written by ChargeMigration or densitypy, instead it uses Pulses3D which does not output a csv format
file and it is expected to remain as such. Most ofther files are csv.

Reading of grid densities for each state was fixed to properly match Molcas dynamic chunk sizes. the entire pipeline was
streamlined  