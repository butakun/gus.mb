// $Id: DataExchanger.h 21 2010-05-17 13:50:02Z kato $
#ifndef INCLUDED_DATA_EXCHANGER_H__
#define INCLUDED_DATA_EXCHANGER_H__

class DataExchanger
{
public:
    virtual ~DataExchanger() {}

    virtual void Start() = 0;
    virtual void Finish() = 0;

protected:

private:
};

#endif // INCLUDED_DATA_EXCHANGER_H__

