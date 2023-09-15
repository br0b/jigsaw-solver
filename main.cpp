#include<algorithm>
#include<cassert>
#include<cstdio>
#include<cstdlib>
#include<vector>

// element of a matrix
// i - num of row
// j - num of column
struct Element {
  int i;
  int j;
};

class Matrix {
  protected:
    // holds elements of matrix
    std::vector<std::vector<char>> matrix;
    // returns first free element (counting from the left)
    // from the first non-empty row (counting from the top)
    Element getFreeElement() {
      for (int i = 0; i < (int)getNumOfRows(); i++)
        for (int j = 0; j < (int)getNumOfColumns(); j++)
          if (get(i, j) == '.') return { i, j };
      exit(1);
    }
    bool rowEmpty(int i) {
      for (char c : matrix[i])
        if (c != '.') return false;
      return true;
    }
    bool rowEmpty(std::vector<char> &row) {
      for (char el : row)
        if (el != '.') return false;
      return true;
    }
    bool columnEmpty(int j) {
      for (auto &row : matrix)
        if (row[j] != '.') return false;
      return true;
    }
  public:
    Matrix() : matrix() {}
    // n - num of rows
    // m - num of columns
    Matrix(int n, int m) : matrix(n, std::vector<char>(m)) {
      for (auto &row : matrix)
        for (char &el : row)
          el = '.';
    }
    bool empty() {
      for (auto &row : matrix)
        if (!rowEmpty(row)) return false;
      return true;
    }
    bool full() {
      for (auto &row : matrix)
        for (char &el : row)
          if (el == '.')
            return false;
      return true;
    }
    int getNumOfRows() const { return (int)matrix.size(); }
    int getNumOfColumns() const { return (int)matrix.back().size(); }
    // return num of elements with symbol different from '.'
    int getNumOfElems() const {
      int elems = 0;
      for (const auto &row : matrix)
        for (const char &c : row)
          if (c != '.') elems++;
      return elems;
    }
    int size() { return getNumOfRows() * getNumOfColumns(); }
    char get(int i, int j) const { return matrix[i][j]; }
    void set(int i, int j, char c) { matrix[i][j] = c; }
    // eliminate the free rows and columns
    void crop() {
      if (empty()) return;
      // upper left corner of the section about to be cut out
      Element ul= { 0, 0 };
      while (rowEmpty(ul.i)) ul.i++;
      while (columnEmpty(ul.j)) ul.j++;
      // lower right corner of the section about to be cut out
      Element lr = ul;
      while (lr.i != getNumOfRows() && !rowEmpty(lr.i)) lr.i++;
      while (lr.j != getNumOfColumns() && !columnEmpty(lr.j)) lr.j++;
      int num_of_rows = lr.i - ul.i;
      int num_of_columns = lr.j - ul.j;
      std::vector<std::vector<char>> cropped(num_of_rows,
        std::vector<char>(num_of_columns));
      for (int i = 0; i < num_of_rows; i++)
        for (int j = 0; j < num_of_columns; j++)
          cropped[i][j] = matrix[ul.i + i][ul.j + j];
      matrix = cropped;
    }
    // rotate 90 degrees clockwise
    void rotate() {
      int new_n_of_rows = (int)getNumOfColumns();
      int new_n_of_cols = (int)getNumOfRows();
      std::vector<std::vector<char>> new_matrix(new_n_of_rows,
        std::vector<char>(new_n_of_cols));
      for (int i = 0; i < new_n_of_rows; i++)
        for (int j = 0; j < new_n_of_cols; j++)
          new_matrix[i][j] = matrix[new_n_of_cols - j - 1][i];
      matrix = new_matrix;
    }
    // set element symbols from standard input
    void setFromStdIn() {
      for (int i = 0; i < getNumOfRows(); i++) {
        for (int j = 0; j < getNumOfColumns(); j++) {
          set(i, j, (char)getchar());
        }
        getchar();
      }
    }
    bool fits(Matrix &other) {
      return getNumOfColumns() <= other.getNumOfColumns() &&
             getNumOfRows() <= other.getNumOfRows();
    }
    void print() {
      for (std::vector<char> row : matrix) {
        for (char el : row)
          putchar(el);
        putchar('\n');
      }
    }
    bool operator==(const Matrix &other) { return matrix == other.matrix; }
    bool operator<(const Matrix &other) const {
      return getNumOfElems() < other.getNumOfElems();
    }
};

class Piece {
  private:
    // different rotations of the same piece
    std::vector<Matrix> rotations = {};
    int rotation = 0;
    // if piece is prioritised, it means that it
    //  - is put on the board as first
    //  - if it is taken out, than it isn't possible
    //    to complete the board with the given set
    bool prioritised = false;
    bool used = false;
  public:
    void set(Matrix matrix) {
      Matrix cropped = matrix;
      cropped.crop();
      // check if matrix fits perfectly to board
      if (cropped.getNumOfColumns() == matrix.getNumOfColumns() &&
          cropped.getNumOfRows() == matrix.getNumOfRows()) {
        // a 180 rotation of perfectly fitting peice is not interesting
        rotations = { cropped };
        prioritised = true;
        return;
      }
      Matrix _rotations[4];
      for (int i = 0; i < 4; i++) {
        _rotations[i] = cropped;
        cropped.rotate();
      }
      //check if matrix is symmetrical
      if (_rotations[0] == _rotations[1])
        rotations = { _rotations[0] };
      else if (_rotations[0] == _rotations[2]) 
        rotations = { _rotations[0], _rotations[1] };
      else
        rotations = { _rotations[0], _rotations[1],
                      _rotations[2], _rotations[3] };
    }
    const Matrix* getMatrix() const { 
      assert(!rotations.empty());
      return  &rotations[0];
    }
    const Matrix* getMatrix(int r) const {
      assert(!rotations.empty());
      return &rotations[r];
    }
    int getNumOfUniqueRotations() { return (int)rotations.size(); }
    int getNumOfElems() {
      assert(!rotations.empty());
      return rotations[0].getNumOfElems();
    }
    // align matrix's upper left corner so it covers element el
    Element align(Element el, int r) {
      int j = 0;
      while (rotations[r].get(0, j) == '.') j++;
      return { el.i, el.j - j };
    }
    bool onBoard() { return used; }
    bool isPrioritised() const { return prioritised; }
    void use() {
      assert(used == false);
      used = true;
    }
    void free() {
      assert(used == true);
      used = false;
    }
    void print() { rotations[rotation].print(); }
    bool operator<(const Piece &other) {
      if (isPrioritised()) return false;
      if (other.isPrioritised()) return true;
      return *(this->getMatrix()) < *(other.getMatrix());
    }
};

// set of pieces
class Set {
  private:
    // this is a subset of the vector this pointer points to
    std::vector<Piece> *pieces;
    // a *(pieces)[i] belongs to this if and only if 
    // i belongs to vector <set>.
    // vector set is sorted in an ascending order
    std::vector<int> set;
    int sumOfSizes() {
      if (set.empty()) return 0;
      assert(pieces != nullptr);
      int sum = 0;
      for (int piece : set) sum += (*pieces)[piece].getNumOfElems();
      return sum;
    }
  public:
    Set() : pieces{nullptr} {}
    Set(std::vector<Piece> &p) : pieces{&p} {}
    // return true if there are no pieces in set, false otherwise
    bool empty() { return set.empty(); }
    // return biggest index currently in set.
    // if empty(), than this causes and undefined behaviour
    int back() {
      return set.back();
    }
    void push(int i) {
      assert(set.empty() || set.back() < i);
      set.push_back(i);
    }
    // return vector of indexes of pieces belonging to this 
    std::vector<int> get() { return set; }
    bool isInteresting(int matrix_size) { 
      return this->sumOfSizes() == matrix_size;
    }
};

// represents an addition one can make to a board
class Addition {
  private:
    int piece_id;
    // upper left corner of piece
    Element upper_left_corner;
    int rotation;
  public:
    Addition(int i, Element el, int r) : piece_id(i),
      upper_left_corner(el),
      rotation(r) {}
    int getPieceId() { return piece_id; }
    Element getUpperLeft() { return upper_left_corner; }
    int getNumOfRotations() { return rotation; }
    bool isInteresting(Matrix *board, std::vector<Piece> &pieces) {
      Matrix const *m = pieces[piece_id].getMatrix(rotation);
      assert(upper_left_corner.i >= 0);
      // check if piece would cross borders of the board
      if (upper_left_corner.j < 0 ||
          upper_left_corner.i + m->getNumOfRows() >
            board->getNumOfRows() ||
          upper_left_corner.j + m->getNumOfColumns() >
            board->getNumOfColumns())
        return false;
      // check if the section of the board to be covered
      // is now free
      for (int i = 0; i < (int)m->getNumOfRows(); i++)
        for (int j = 0; j < (int)m->getNumOfColumns(); j++)
          if (board->get(upper_left_corner.i + i,
                         upper_left_corner.j + j) != '.' &&
              m->get(i, j) != '.')
            return false;
      return true;
    }
};

class Board : public Matrix {
  private:
    // return a vector of subsets of pieces that are interesting, as in whose
    // sum of non-empty elements is equal to the size of the board
    std::vector<Set> interestingSubsets(std::vector<Piece> &pieces) {
      const Set EMPTY(pieces);
      // group subsets with the same number of pieces together
      std::vector<std::vector<Set>> subsets(pieces.size() + 1);
      // vector of interesting subsets to return
      std::vector<Set> interesting_subsets;
      // start constructing pieces from an empty piece
      subsets[0] = { EMPTY };
      for (int num_of_el = 1; num_of_el <= (int)subsets.size(); num_of_el++)
        for (Set &subset : subsets[num_of_el - 1]) {
          int starting_piece;
          if (subset.empty()) starting_piece = 0;
          else starting_piece = subset.back() + 1;
          // create a new subset by adding a single piece to an already
          // created one. The new index is bigger than the biggest in
          // the already created set
          for (int piece = starting_piece; piece < (int)pieces.size();
               piece++) {
            Set new_subset = subset;
            new_subset.push(piece);
            subsets[num_of_el].push_back(new_subset);
            if (new_subset.isInteresting(this->size())) {
              interesting_subsets.push_back(new_subset);
            }
          }
        }
      return interesting_subsets;
    }
    std::vector<Addition> possibleAdditions(Set &subset,
      std::vector<Piece> &pieces) {
      Element free_el = getFreeElement();
      std::vector<Addition> additions;
      for (int i : subset.get()) {
        if (pieces[i].onBoard()) continue;
        // put the prioritised board in a random way, but legal way.
        // it is the only possible move
        if (pieces[i].isPrioritised()) return { Addition(i, { 0, 0 }, 0) };
        for (int r = 0; r < pieces[i].getNumOfUniqueRotations(); r++) {
          // align upper left corner so that the free element is covered
          Element ul = pieces[i].align(free_el, r);
          Addition addition(i, ul, r);
          if (addition.isInteresting(this, pieces))
            additions.push_back(addition);
        }
      }
      return additions;
    }
    void makeAddition(Addition &add, std::vector<Piece> &pieces) {
      assert(pieces[add.getPieceId()].onBoard() == false);
      pieces[add.getPieceId()].use();
      const Matrix *m =
        pieces[add.getPieceId()].getMatrix(add.getNumOfRotations());
      for (int i = 0; i < (int)m->getNumOfRows(); i++)
        for (int j = 0; j < (int)m->getNumOfColumns(); j++) {
          int board_i = add.getUpperLeft().i + i;
          int board_j = add.getUpperLeft().j + j;
          assert(0 <= board_i && board_i < getNumOfRows());
          assert(0 <= board_j && board_j < getNumOfColumns());
          assert(get(board_i, board_j) == '.' ||
                 m->get(i, j) == '.');
          if (m->get(i, j) != '.') set(board_i, board_j, m->get(i, j));
        }
    }
    void revertAddition(Addition &add, std::vector<Piece> &pieces) {
      assert(pieces[add.getPieceId()].onBoard());
      pieces[add.getPieceId()].free();
      const Matrix *m =
        pieces[add.getPieceId()].getMatrix(add.getNumOfRotations());
      for (int i = 0; i < (int)m->getNumOfRows(); i++)
        for (int j = 0; j < (int)m->getNumOfColumns(); j++) {
          int board_i = add.getUpperLeft().i + i;
          int board_j = add.getUpperLeft().j + j;
          assert(0 <= board_i && board_i < getNumOfRows());
          assert(0 <= board_j && board_j < getNumOfColumns());
          if (m->get(i, j) != '.') {
            assert(get(board_i, board_j) == m->get(i, j));
            set(board_i, board_j, '.');
          }
        }
    }
    void _fill(Set &subset, std::vector<Piece> &pieces) {
      if (full()) return;
      for (Addition add : possibleAdditions(subset, pieces)) {
        makeAddition(add, pieces);
        _fill(subset, pieces);
        if (full()) return;
        revertAddition(add, pieces);
      }
    }
  public:
    Board(int a, int b) : Matrix(a, b) {}
    void fill(std::vector<Piece> pieces) {
      // sort the pieces so the bigger ones are upfront.
      // They will be considered first
      std::sort(pieces.begin(), pieces.end(),
        [](Piece a, Piece b){ return b < a; });
      for (Set &subset : interestingSubsets(pieces)) {
        _fill(subset, pieces);
        if (full()) return;
        assert(empty());
      }
    }
};

int main() {
  int n;
  int m;
  int num_of_pieces;
  assert(scanf("%d %d %d", &n, &m, &num_of_pieces));
  getchar();
  if (num_of_pieces == 0) {
    printf("NIE");
    return 0;
  }
  Board board(n, m);
  std::vector<Piece> pieces(num_of_pieces);
  Matrix user_input(n, m);
  user_input.setFromStdIn();
  pieces[0].set(user_input);
  for (int i = 1; i < num_of_pieces; i++) {
    // read the additional line between brick descriptions
    assert(getchar() == '\n');
    user_input.setFromStdIn();
    pieces[i].set(user_input);
  }
  board.fill(pieces);
  if (board.full()) {
    printf("TAK\n");
    board.print();
    return 0;
  }
  printf("NIE");
}
