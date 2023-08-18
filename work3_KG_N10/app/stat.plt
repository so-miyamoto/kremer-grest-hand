
set terminal png
set output "result.png"
set key right outside
f="stat.log"
p f u 1:2 t "potential Eall"   ,\
  f u 1:3 t "potential Ebond"  ,\
  f u 1:4 t "potential Epair"  ,\
  f u 1:5 t "kinetic Ene K"    ,\
  f u 1:6 t "all Ene"          ,\
  f u 1:7 t "temperature"      ,\
  f u 1:8 t "pressure"    

