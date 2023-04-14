version 1.0

struct Person {
  String name
  Int age
}

task apps_1381 {
  input {
    Person person
    Int years
  }

  command {
    echo "hello my name is ${person.name} and I am ${person.age} years old, although writing WDL has aged me ~{years} years"
  }

  output {
    Person older_person = object {
      name: person.name,
      age: person.age + years
    }
    String message = read_string(stdout())
  }
}