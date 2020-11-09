#ifndef TIME_SPENT_H
#define TIME_SPENT_H

//#include <sys/time.h>
#include <string>
#include <ctime>
#include <ratio>
#include <chrono>

class TimeSpent
{
public:
    //! construtor
    TimeSpent(std::string = "");

    //! destrutor
    ~TimeSpent();

    //! sets
    void setName(std::string _n) {name = _n;}

    //! começa a contagem do tempo
    void startCount();

    //! finaliza a contagem do tempo
    void endCount();

    //! imprime o tempo
    void print();

    void setZero();

    std::chrono::duration<double> getAccumulator() const;

private:
    std::string name;
    std::chrono::system_clock::time_point start, end;
    std::chrono::duration<double> accumulator;
};

#endif // TIME_SPENT_H
