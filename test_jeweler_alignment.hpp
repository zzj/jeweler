#include <string>
#include <memory>
using namespace std;
class JewelerAlignment;

void create_alignment(const int position, const string &cigar_string,
                      const char rep, int query_length, JewelerAlignment &al);

class CreateAlignment {
public:
    CreateAlignment();
    CreateAlignment& set_position(int position);
    CreateAlignment& set_length(int length);
    CreateAlignment& set_query_bases(char rep);
    CreateAlignment& set_qualities(char rep);
    CreateAlignment& set_cigarop(string cigar_string);
    JewelerAlignment get();

private:
    unique_ptr<JewelerAlignment> al;
    int length;
};
