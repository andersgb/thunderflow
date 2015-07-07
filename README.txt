Install instructions:
- depends on rgrad, clone rgrad into cpp/rgrad
- put pr_c4.rb cpp/rgrad

- rgrad requires ruby
- requires boost-python

cd cpp/rgrad
ruby rgrad.rb --generate=C --outputdir=models pr_c4.rb
cd ../..
make -C cpp
python run.py
python view.py 
