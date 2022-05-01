#!/usr/bin/env sage
from random import choice, sample
import json
import signal
import datetime

r"""
    An example of how to use the class "Sbox" from external code

AUTHORS:

- Oleksandr Kazymyrov (2011-06-21): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013 Oleksandr Kazymyrov <okazymyrov@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

os.chdir(os.path.split(os.path.abspath(sys.argv[0]))[0] + "/Sage")

load("./TestFunctions.sage")

def gen_new_sbox(method='random_substitution'):
    S=Sbox(n=8,m=12)
    S.generate_sbox(method=method)
    return S

def permute_sbox(sbox):    
    out_sbox = Sbox(n=8,m=12)
    idx = range(len(sbox))
    i1, i2 = sample(idx, 2)
    sbox[i1], sbox[i2] = sbox[i2], sbox[i1]
    out_sbox.set_sbox(sbox)
    return out_sbox


def assess_sbox_solution(sbox):
    ret=sbox.algebraic_immunity_sbox()
    if sbox.fixed_points():
        score = 0
    else:
        score = (sbox.NL() * 0.01) \
            + (sbox.minimum_degree() * 0.15) \
            + ((ret[0]+(ret[1]/1000))*0.33) \
            + (sbox.MDT()*-0.125)
    return {"sbox": sbox.get_sbox(), "minimum_degree":sbox.minimum_degree(), "algebraic_immunity": ret[0], "number_algebraic_equations": ret[1], "nonlinearity": sbox.NL(), "uniformity": sbox.MDT(), "fixed_points": sbox.fixed_points(), "score": score }

def generate_next_sbox(solutions):
    
    # If we have existing solutions we can choose to modify an existing solution from the top 10
    solutions.sort(key=lambda sol:sol['score'], reverse=True)
    start_sbox = choice(solutions[:10]) if solutions else None
    if start_sbox:
        # try 100 permutations to find a better solution
        for i in range(100):
            sbox = permute_sbox(start_sbox['sbox'])
            sbox = assess_sbox_solution(sbox)
            if sbox['score'] > start_sbox['score']:
                print(f"Permuted better solution in {i} iterations")
                return sbox
            start_sbox = sbox

    # return a newly generated sbox
    print("Generating new random solution")
    return assess_sbox_solution(gen_new_sbox())

def save_results(solutions):
    with open(f"results/{datetime.datetime.utcnow().isoformat().replace(':','-')}--{len(solutions)}.json","w") as f:
        json.dump(solutions, f, cls=json_sage)
        
class json_sage(json.JSONEncoder):
    def default(self, o):
        try:
            if type(o)==sage.rings.integer.Integer:
                return int(o)
            elif type(o)==sage.rings.real_mpfr.RealLiteral:
                return float(o)
            elif type(o)==sage.rings.real_mpfr.RealNumber:
                return float(o)
            elif type(o)==sage.symbolic.expression.Expression:
                return latex(o)
        except TypeError:
            pass

        # Let the base class default method raise the TypeError
        return json.JSONEncoder.default(self, o)

def signal_handler(signal, frame):
    global interrupted
    interrupted = True

signal.signal(signal.SIGINT, signal_handler)

interrupted = False

def main(argv=None):
    global interrupted
    # bits=8
    solutions = []
    solution_count = 0

    t1=cputime()
    t2=cputime()
    # test_all_functions(bits=bits)
    # test_temp(bits=bits)
    # gen_8_12()
    print("Seeding solutions")
    for generator_method in ['dickson','gold','inverse','kasami','niho','random_substitution','welch']:    
        solution = assess_sbox_solution(gen_new_sbox(generator_method))
        print(f"Seeding solution with score {solution['score']}, minimum_degree {solution['minimum_degree']}, algebraic_immunity {solution['algebraic_immunity']}, nonlinearity {solution['nonlinearity']}, uniformity {solution['uniformity']}" )
        solutions.append(solution)
    
    # Keep looking for solutions until a specific amount of time has passed or until we have 8    
    while ((t2 - t1) < 60*60*5 or len(solutions) < 8) and not interrupted:  
        solution_count+=1
        solution = generate_next_sbox(solutions)
        print(f"Found solution with score {solution['score']}, minimum_degree {solution['minimum_degree']}, algebraic_immunity {solution['algebraic_immunity']}, nonlinearity {solution['nonlinearity']}, uniformity {solution['uniformity']}" )
        solutions.append(solution)        
        save_results(solutions)
        t2=cputime()

    print(f"Found {len(solutions)} solutions.")

    print("=====")
    print("Time = {0}".format(t2-t1))
    print("Seconds per solution = {0}".format((t2-t1)/solution_count))


if __name__ == "__main__":
    sys.exit(main())
