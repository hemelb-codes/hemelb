set breakpoint pending on

# break Neighbours::Init()
break ConfigGenerator::Execute()

# break LaterNeighbourIterator::LaterNeighbourIterator(Site& site)
# break LaterNeighbourIterator::LaterNeighbourIterator(Site& site, unsigned int startpos)
# break LaterNeighbourIterator::LaterNeighbourIterator(const LaterNeighbourIterator& other)
# break LaterNeighbourIterator::GetVector()
# break LaterNeighbourIterator::IsCurrentValid()
# break LaterNeighbourIterator::operator++()
# break LaterNeighbourIterator::operator*()

# break ConfigGenerator.cpp:215
break BlockWriter::Finish()
