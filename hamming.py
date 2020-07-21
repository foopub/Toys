import random
from math import *


def hamming(string: str):
    """
    Return a set of bad parity bits for a given string.
    """
    bad = set()

    for i in range(floor(log(len(string), 2))+1):
        x = ''
        start = 2**i-1
        end = 2**(i+1)-1

        while start < len(string):
            try:
                x += string[start:end]
            except:
                pass
            start += 2**(i+1)
            end += 2**(i+1)

        if x[1:].count('1') % 2 != int(string[2**i-1]):
            bad.add(i)  

    return bad


def errorbit(bad: set):
    """
    Return the position of an erroneous bit, from set of bad parities.
    """
    error = -1

    for i in bad:
        error += 2**i
    
    return error


def parities(safe: int):
    """
    Express a number as a sum of powers of two.
    """
    result = set()

    while safe:
        n = floor(log(safe, 2))
        result.add(n)
        safe -= 2**n

    return result


def prisoner(string: str):
    """
    Return the position of an erroneous bit in a string.
    """

    return errorbit(hamming(string))+1
    
    
def change_bits(string: str, *bits: int):
    """
    Flip the bits at the positions specified.
    """

    for i in bits:
        bit = str(int(int(string[i]) + 1 == 1))
        string = string[:i] + bit + string[i+1:]

    return string


def guard(string: str, safe: int):
    """
    Return a string with error in specified position, by changing <= 2 bits.
    """
    old = hamming(string)
    error = errorbit(old)
    new = parities(safe)
    change = new ^ old
    bit = errorbit(change)

    print(f'You know the switch is #{safe} of', len(string)-1)
    print(f'The current state is:\n{old}')
    print(f'Desired state is:\n{new}')
    print(f'Therefore change:\n{change}')

    if not change:
        print('...nothing to change')
        return string
    elif bit < len(string):
        print(f'by changing bit:\n{bit}')
        return change_bits(string, bit)
    else:
        bit1 = 2**floor(log(bit, 2))-1
        bit2 = bit - bit1-1
        print(f'But because {bit} is outside the range, change {bit1} and {bit2}')
        return change_bits(string, bit1, bit2)


def generate(N):
    #switch is between 0 and N
    correct_switch = random.randint(0, N-1)
    #one bit for every switch
    bitlist = []
    for i in range(0,N):
        x = str(random.randint(0, 1))
        bitlist.append(x)

    sequence = ''.join(bitlist)

    return (correct_switch, sequence)

def differences(sequence1, sequence2):
    differences = 0
    for i in range(0, len(sequence1)):
        if sequence1[i] is not sequence2[i]:
            differences += 1
    print(differences)
    return differences

def main():
    (switch, sequence) = generate(random.randint(2, 100))
    print("Switch:%d, sequence:%s" % (switch, sequence))
    new_sequence = guard(sequence, switch)
    print("Guard sequence:%s" % (new_sequence))
    guessed_switch = prisoner(new_sequence)                             #small mistake here
    print("Prisoner guess:%d" % guessed_switch)
    #make sure prisoner picks right switch
    assert switch == guessed_switch
    #make sure no extra switches are added in the process
    assert len(new_sequence) <= len(sequence)
    #make sure no more than 4 switches are being flipped
    assert differences(sequence, new_sequence) <= 2   #changed to 2

count = 0
while count < 1000:
        main()
        count += 1
        print(f'Current count is {count}')
